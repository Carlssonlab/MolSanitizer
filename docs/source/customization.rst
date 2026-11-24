.. _customization:

Customization
=============

MolSanitizer reads its chemistry rules from plain-text files that ship with the
package. Each file can be replaced by your own version:

.. list-table::
   :header-rows: 1
   :widths: 22 22 26 30

   * - Rules
     - Use your file
     - Write the template
     - Template file name
   * - Salts and fragments
     - ``--saltlib``
     - ``--create_saltlib``
     - ``salt_stripping.txt``
   * - Protonation
     - ``--protlib``
     - ``--create_protlib``
     - ``ionizations_v3.txt``
   * - Tautomers
     - ``--taulib``
     - ``--create_taulib``
     - ``tautomers_v3.txt``
   * - Substructure filters
     - ``--custom``
     - ``--create_custom``
     - ``filter_out.txt``

The ``--create_*`` options copy the bundled file into the current working
directory (with ``--prefix NAME``, the copy is written as ``NAME.txt``).

A file passed with ``--saltlib``, ``--protlib`` or ``--taulib`` **replaces** the
bundled file completely. To extend the defaults rather than replace them, start
from the template and keep every rule that should remain active.

``--custom`` works differently: it does not replace anything. Its patterns are
applied as an additional filter, together with the built-in filters selected by
``--unwanted`` if that option is also given.

In all files, lines starting with ``#`` and empty lines are ignored. Fields are
separated by whitespace (tabs or spaces), so a field must not contain spaces
itself.

Protonation rules
-----------------

Each line holds six fields:

.. code-block:: text

    FUNCTIONAL_GROUP   pKa   TYPE   Enumerate   REACTION   REF

Only the first five are required.

* ``FUNCTIONAL_GROUP`` is the rule name (without spaces), usually the functional group being ionized.

* ``pKa`` is the estimated pKa value of the functional group, given as an integer (without decimals).
  Only one value per group can be assigned.

* ``TYPE`` is ``ACID`` (deprotonation) or ``BASE`` (protonation).
* An ``ACID`` rule is applied when its pKa is at or below the pH, and a
  ``BASE`` rule when its pKa is at or above the pH. Rules are evaluated at the
  chosen ``--pH`` and across the range set by ``--pH_range``.

* ``Enumerate = 0``: the rule is applied at one matching site, then repeated on
  the product until no site matches, so all matching sites react.
* ``Enumerate = 1``: each matching site starts a separate branch. If the SMARTS
  prevents a second site from reacting (for example, by excluding sites close
  to an existing charge, as the default amine rule does), each site yields its
  own product. For example, 1-methylpiperazine at pH 7 gives two monocations
  and no dication:

  .. code-block:: console

      $ msani -s "CN1CCNCC1" --protonation --pH 7
      C[NH+]1CCNCC1 1_1
      CN1CC[NH2+]CC1 1_2

* ``REACTION`` is a reaction SMARTS.

The rules run in this order:

1. First all ``ACID`` rules by ascending pKa,
   then all ``BASE`` rules by descending pKa.
2. Rules with the same pKa: first ``Enumerate = 0``,
   then ``Enumerate = 1``.
3. Rules that are still tied: in file order.

Tautomer rules
--------------

Each line holds four fields:

.. code-block:: text

    name    enumerate	smarts      reference

Only the first three are required.

* ``name`` (without spaces), usually the functional group being tautomerized.

* ``enumerate = 0`` **standardizes**: the input is replaced by the preferred
  tautomer. For example, an enol is converted to the ketone
  (``CC(O)=CC`` → ``CCC(C)=O``).
* ``enumerate = 1`` **enumerates**: the input is kept and alternative forms are
  added. For example, both the 1H and 3H tautomers of 4-methylimidazole are
  returned.

* ``smarts`` is the reaction pattern in SMARTS format.

The rules run in this order:

1. By default, the input is first canonicalized with RDKit's tautomer
   canonicalizer (``--notaurdkit`` turns this off).
2. All ``Enumerate = 0`` rules, in file order.
3. All ``Enumerate = 1`` rules, in file order.

``--extended-tautomers`` adds a second bundled file (``tautomers_extended.txt``)
with rules for less probable forms. These rules run after your enumerating
rules and are **not** replaced by ``--taulib``. ``--extended-tautomers`` also
turns on ``--tautomers``.

Salts and fragments
-------------------

Each line holds a SMARTS pattern, optionally followed by a name:

.. code-block:: text

    [Cl] chloride
    [CH3]C(=O)O acetate

Fragments matching a pattern are removed. If more than one fragment remains
afterwards, whether or not any pattern matched, only the largest fragment (by
heavy-atom count) is kept. If every fragment matches a pattern, the entry is
removed.

``--saltlib`` requires ``--removesalts`` (or ``--lazy``, which turns it on).

For example, with the acetate line removed from the salt list, sodium acetate
(``CC(=O)O.[Na]``) keeps its acetic acid: only the sodium is stripped. With
the default list, both fragments match and the entry is removed.

Substructure filters
--------------------

The first two columns, ``SMARTS`` and ``LABEL``, are required; further columns
are ignored. A molecule matching any pattern is rejected and written to the
rejected-entries file with its label. See :doc:`usage` for the built-in ``--unwanted`` categories.

Worked example: a new protonation rule
--------------------------------------

As an illustration, we add a rule that deprotonates aliphatic alcohols at pH 7.
This is chemically unrealistic and only shows how a rule works. The rule line
has five fields:

.. code-block:: text

    alcohols	5	ACID	0	[OH1&+0;$(O-[CX4&z1])!$(O~*~[*-]):1]>>[OH0&-:1]

1. ``alcohols``: the rule name, without spaces.
2. ``5``: the estimated pKa.
3. ``ACID``: the type, ``ACID`` for deprotonation and ``BASE`` for protonation.
4. ``0``: the ``Enumerate`` flag, so the rule is applied at all matching sites.
5. The reaction SMARTS. It selects a neutral hydroxyl oxygen bonded to an sp3
   carbon (``CX4``, four connections in total) that has no other heteroatom
   neighbour (``z1``). ``!$(O~*~[*-])`` skips a hydroxyl whose neighbouring
   atom is already bonded to a negatively charged atom.

The SMARTS reuses patterns from the bundled file, for example
``[OH1&+0;...`` and ``...>>[OH0&-:1]`` from the OH-acid rules.

The line can be saved on its own in a new file, in which case only this rule
is active, or appended to a copy of the template written by
``--create_protlib``, in which case it is used in addition to the default
rules. Save it as ``aliphatic_alcohol.txt`` and run:

.. code-block:: console

    $ msani -s "CO" --protonation --protlib aliphatic_alcohol.txt
    C[O-] 1

Methanol and tert-butanol (``CC(C)(C)O``) are converted to the alkoxides
``C[O-]`` and ``CC(C)(C)[O-]``. With the default library, both remain neutral.

Usage
-----

.. code-block:: console

    $ msani --create_protlib
    $ msani -i input.smi --protonation --pH 7 --protlib my_prot.txt --tautomers --taulib my_tau.txt
    $ msani -i input.smi --removesalts --saltlib my_salts.txt

The same options can be set in a YAML file passed with ``--config``
(see :ref:`config`):

.. code-block:: yaml

    protlib: my_prot.txt
    taulib: my_tau.txt
    saltlib: my_salts.txt
    custom: my_filters.txt

The standalone tools ``m_ionizer`` and ``m_tautomerizer`` take a rules file
with ``--library`` (``-l``).

From Python, pass ``salt_library=``, ``protonation_library=`` and
``tautomer_library=`` to :class:`~msani.api.Msani`, or
``salt_library=`` to :class:`~msani.filtering.filters.Filters`.

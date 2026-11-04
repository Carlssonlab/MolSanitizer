/* Keep the release number visible alongside Read the Docs' version alias. */
(() => {
  function showReleaseVersion() {
    const label = document.querySelector(".msani-doc-version");
    const selector = document.querySelector(".version-switch select");
    if (!label || !selector) {
      return;
    }

    selector.setAttribute("aria-label", "Documentation version");
    const option = selector.selectedOptions[0];
    if (option && (option.value === "stable" || option.value === "latest")) {
      // Keep the selector compact; the full commit suffix is in the label below.
      const publicVersion = label.dataset.release.split("+")[0];
      option.textContent = `${publicVersion} (${option.value})`;
      selector.title = `${label.dataset.release} (${option.value})`;
    }
  }

  document.addEventListener("readthedocs-addons-data-ready", () => {
    // The theme creates the selector during this event. Update it afterward,
    // regardless of which script registered its listener first.
    queueMicrotask(showReleaseVersion);
  });
  document.addEventListener("DOMContentLoaded", showReleaseVersion);
  showReleaseVersion();
})();

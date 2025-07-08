from .geometry import distL2Squared3
import itertools
import math

class buckets2:

		def __init__(self, tolerance):

			self.tolerance = tolerance
			self.tolerance2 = tolerance**2
			# see explanation for bucketsize below
			self.bucketsize = tolerance / math.sqrt(3)
			self.extrawidth = tolerance - self.bucketsize

			#print(self.bucketsize, self.extrawidth)

			self.moore_neighborhood = list(itertools.product(range(-1, 2), range(-1, 2), range(-1, 2)))
			self.moore_neighborhood.remove((0, 0, 0))

			self.moore_faces = [None for i in range(6)]
			self.moore_faces[0] = [(-2, i, j) for i,j in itertools.product(range(-1, 2), range(-1, 2))]
			self.moore_faces[1] = [( 2, i, j) for i,j in itertools.product(range(-1, 2), range(-1, 2))]
			self.moore_faces[2] = [(i, -2, j) for i,j in itertools.product(range(-1, 2), range(-1, 2))]
			self.moore_faces[3] = [(i,  2, j) for i,j in itertools.product(range(-1, 2), range(-1, 2))]
			self.moore_faces[4] = [(i, j, -2) for i,j in itertools.product(range(-1, 2), range(-1, 2))]
			self.moore_faces[5] = [(i, j,  2) for i,j in itertools.product(range(-1, 2), range(-1, 2))]

		def _get_moore_faces(self, xyz, bxyz):

			x, y, z = xyz
			bx, by, bz = bxyz
			xd, yd, zd = (x-(bx*self.bucketsize)), (y-(by*self.bucketsize)), (z-(bz*self.bucketsize)) 

			faces = []

			if xd < self.extrawidth:
				faces.append(0)
			if xd > (self.bucketsize - self.extrawidth):
				faces.append(1)
			if yd < self.extrawidth:
				faces.append(2)
			if yd > (self.bucketsize - self.extrawidth):
				faces.append(3)
			if zd < self.extrawidth:
				faces.append(4)
			if zd > (self.bucketsize - self.extrawidth):
				faces.append(5)

			return faces

		def bucket(self, xyzData, atomId=0, confNum=0, confClusters=None):

			natoms = len(xyzData)
			buckets = {}
			npos = 0

			for atom in range(natoms):
				
				xyzactual = xyzData[atom]
				bx, by, bz = math.floor(xyzactual[0]/self.bucketsize), math.floor(xyzactual[1]/self.bucketsize), math.floor(xyzactual[2]/self.bucketsize)
				xyzhash = (bx, by, bz)
				visited = [False]
				if not buckets.get(xyzhash):
					buckets[xyzhash] = (visited, [atom]) # using a list so I can treat the visited flag as a pointer, since tuples are immutable
				else:
					buckets[xyzhash][1].append(atom)

			for xyzhash, (visited_this, atoms) in buckets.items():
				if not atoms: # if all mols from this cluster have already been absorbed into another cluster we continue
					continue

				xyz = xyzData[atoms[0]]

				currCluster = [m for m in atoms] 
				# any mols in this bucket are automatically clustered together, since by definition they must be within the tolerance distance
				# bucket_size := (tolerance) / sqrt(3)
				# ( this is the maximal size of a cube inscribed in a sphere where diameter == tolerance ) 
				# Since the diagonal length of a cube is sqrt(3) * r, even if two points are spaced out maximally within the cube at opposite ends of a diagonal, they will be within the tolerance distance
				# problem is that this requires checking an extended moore neighborhood of neighboring buckets. From 26 to 80 neighbors.
				# however, by analyzing the exact position of the point in each bucket, we can reduce the number of neighbors to check from an average of 80 to an average of ~65

				# may want to scrap the idea of buckets being guaranteed clusters in the future, because it requires more hash lookups
				# for now this solution is fast enough and robust for getting correct clustering

				def check(loc):
					neighbor_mol2 = buckets.get(loc)

					if neighbor_mol2:
						visited, neighboring_atoms = neighbor_mol2
						if visited[0]:
							return

						lneighbor = len(neighboring_atoms)
						j = 0
						while j < lneighbor:
							ni = neighboring_atoms[j]
							xyzn = xyzData[ni]

							dist = distL2Squared3(xyz, xyzn)
							if dist <= self.tolerance2:
								currCluster.append(ni)
								neighboring_atoms.pop(j)
								j -= 1
								lneighbor -= 1

							j += 1

				# this is where we analyze the point's xyz to determine which additional faces to process
				faces = self._get_moore_faces(xyz, xyzhash)

				# check the immediate moore neighborhood
				for n in range(26):
					neighbor = self.moore_neighborhood[n]
					neighbor_loc = (xyzhash[0]+neighbor[0],xyzhash[1]+neighbor[1],xyzhash[2]+neighbor[2])
					check(neighbor_loc)

				# check the additional faces (at least 3, if the point is in one of the cube corner areas)
				for face_index in faces:
					for n in range(9):
						neighbor = self.moore_faces[face_index][n]
						neighbor_loc = (xyzhash[0]+neighbor[0],xyzhash[1]+neighbor[1],xyzhash[2]+neighbor[2])
						check(neighbor_loc)

				visited_this[0] = True

				# calculate conf clusters as we go
				if confClusters != None:
					currCluster = sorted(currCluster)
					clusterTuple = tuple(currCluster)
					existingcluster = confClusters.get(clusterTuple)
					if not existingcluster:
						confClusters[clusterTuple] = (confNum, [atomId], [xyz])
						confNum += 1
					else:
						confNum_t, atomlist, xyzlist = confClusters[clusterTuple]
						atomlist.append(atomId)
						xyzlist.append(xyz)

				npos += 1

			return npos, confNum
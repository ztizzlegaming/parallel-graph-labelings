# Helper program which takes

import sys

def main():
	# Get the graph parameters and the perm index to graph
	cycleSize = sys.argv[1]
	connectingVertices = sys.argv[2]
	permIdx = int(sys.argv[3])

	# Open the output file
	f = open('output_' + cycleSize + '_' + connectingVertices + '.txt')

	# Convert them to ints
	cycleSize = int(cycleSize)
	connectingVertices = int(connectingVertices)
	
	# First line is the graph parameters, don't need since we got from cmd line
	f.readline()

	# Next line is time taken, not used for anything
	f.readline()

	# Init the adjacency matrix and read it in from output file
	matrix = []
	matrixLine1 = f.readline().rstrip().split()
	vertices = len(matrixLine1)
	matrix.append(matrixLine1)

	for i1 in range(1, vertices):
		matrixLine = f.readline().rstrip().split()
		matrix.append(matrixLine)

	# Convert each edge to an int
	for i1 in range(vertices):
		for i2 in range(vertices):
			matrix[i1][i2] = int(matrix[i1][i2])

	# Get the permutation and convert it to a list
	permLine = None
	for i1 in range(permIdx):
		permLine = f.readline().rstrip()
	permLineParts = permLine.split(': ')
	perm = permLineParts[1][:-13] # Chop off 'Magic Number' off of the end
	perm = perm[1:-1]
	perm = perm.split(', ')

	# Offsets if some vertices need to be higher than others
	outerOffset = 0
	innerOffset = 0
	if connectingVertices > cycleSize - connectingVertices:
		outerOffset = 0.75 * (connectingVertices - (cycleSize - connectingVertices))
	elif connectingVertices < cycleSize - connectingVertices:
		innerOffset = 0.75 * ((cycleSize - connectingVertices) - connectingVertices)

	# Print out the tikz latex code
	print('\\begin{tikzpicture}')
	count = 1

	# Outside vertices on left
	for i1 in reversed(range(cycleSize - connectingVertices)):
		print('\\node[shape=circle,draw=black] (' + str(count) + ') at (0, ' + str(i1 * 1.5 + outerOffset) + ') {' + perm[count - 1] + '};')
		count += 1

	# Inner connecting vertices
	for i1 in range(connectingVertices):
		print('\\node[shape=circle,draw=black] (' + str(count) + ') at (1.5, ' + str(i1 * 1.5 + innerOffset) + ') {' + perm[count - 1] + '};')
		count += 1

	# Outside vertices on right
	for i1 in reversed(range(cycleSize - connectingVertices)):
		print('\\node[shape=circle,draw=black] (' + str(count) + ') at (3, ' + str(i1 * 1.5 + outerOffset) + ') {' + perm[count - 1] + '};')
		count += 1

	# Edges
	for i1 in range(vertices):
		for i2 in range(vertices):
			edge = matrix[i1][i2]
			if edge > 0:
				print('\\path[->] (' + str(i1 + 1) + ') edge node {' + perm[count - 1] + '} (' + str(i2 + 1) + ');')
				count += 1

	print('\\end{tikzpicture}')

	f.close()

if __name__ == '__main__':
	main()
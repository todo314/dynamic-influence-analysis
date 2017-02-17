#include "dim.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	DIM dim;
	dim.init(); // Call at the beginning

	dim.insert(0); // Insert vertex 0
	dim.insert(1); // Insert vertex 1
	dim.insert(0, 1, 0.4); // Insert edge (0,1) with probability 0.5
	dim.set_beta(128); // Set beta=32

	dim.insert(2); // Insert vertex 2
	dim.insert(3); // Insert vertex 3

	dim.insert(1, 2, 0.5); // Insert edge (1,2) with probability 0.6
	dim.insert(2, 0, 0.6); // Insert edge (2,0) with probability 0.7
	dim.insert(2, 3, 0.7); // Insert edge (2,3) with probability 0.8

	for (int v = 0; v < 4; v++)
		printf("Influence of %d is %1.6f\n", v, dim.infest(v));
	printf("Most influential vertex is %d\n", dim.infmax(1)[0]);

	dim.insert(3, 1, 0.8); // Insert edge (3,1) with probability 0.8
	dim.erase(2, 3); // Delete edge (2,3)

	for (int v = 0; v < 4; v++)
		printf("Influence of %d is %1.6f\n", v, dim.infest(v));
	printf("Most influential vertex is %d\n", dim.infmax(1)[0]);

	dim.erase(2); // Delete vertex 2
	dim.change(0, 1, 0.9); // Change probability of edge (0,1) to 0.9

	for (int v = 0; v < 4; v++)
		printf("Influence of %d is %1.6f\n", v, dim.infest(v));
	printf("Most influential vertex is %d\n", dim.infmax(1)[0]);

	return 0;
}

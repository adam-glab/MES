#pragma once

struct Node {
	double x, y;
	bool BC;
	double t0;
	Node() {
		x = 0.; 
		y = 0.; 
		BC = false; 
		t0 = 0.; 
	};
};
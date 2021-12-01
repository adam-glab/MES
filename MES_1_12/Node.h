#pragma once

struct node {
	double x, y;
	bool BC;
	double t0;
	node() { x = 0; y = 0; BC = false; t0 = 0.; };
	node(double x0, double y0) {
		x = x0;
		y = y0;
	}
};
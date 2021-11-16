#pragma once

struct node {
	double x, y;
	bool BC;
	node() { x = 0; y = 0; BC = false; };
	node(double x0, double y0) {
		x = x0;
		y = y0;
	}
};
#pragma once

struct node {
	double x, y;
	node() { x = 0; y = 0; };
	node(double x0, double y0) {
		x = x0;
		y = y0;
	}
};
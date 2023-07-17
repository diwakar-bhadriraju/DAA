#include <iostream>
#include <stack>
#include <algorithm>

struct Point {
    int x, y;
    Point(int a, int b) : x(a), y(b) {}
};

bool compare(Point p1, Point p2) {
    return p1.x < p2.x;
}
int orientation(Point p, Point q, Point r) {
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);
    if (val == 0) {
        return 0;
    }
    return (val > 0) ? 1 : 2;
}

void printConvexHull(Point points[], int n) {
    std::stack<Point> S;
    S.push(points[0]);
    S.push(points[1]);
    for (int i = 2; i < n; i++) {
        while (S.size() > 1 && orientation(S.top(), S.top(), points[i]) != 2) {
            S.pop();
        }
        S.push(points[i]);
    }
    while (!S.empty()) {
        Point p = S.top();
        std::cout << "(" << p.x << ", " << p.y << ")" << std::endl;
        S.pop();
    }
}

int main() {
    const int n = 3;
    Point points[n] = { {1, 2}, {3, 4}, {5, 6} };

    std::sort(points, points + n, compare);

    for (int i = 0; i < n; i++) {
        std::cout << "(" << points[i].x << ", " << points[i].y << ")" << std::endl;
    }

    return 0;
}

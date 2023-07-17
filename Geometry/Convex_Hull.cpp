#include <iostream>
#include <cstdio>
#include <algorithm>
#include <stack>
#include <cmath>

using namespace std;

const double EPS = 1e-9;

class Point {
public:
    double x, y;
    Point() {}
    Point(double x, double y) : x(x), y(y) {}
    bool operator < (const Point& other) const {
        if (fabs(x - other.x) > EPS) {
            return x < other.x;
        } else {
            return y < other.y;
        }
    }
};

int n;
Point points[100005];

double cross(Point a, Point b, Point c) {
    double x1 = b.x - a.x;
    double y1 = b.y - a.y;
    double x2 = c.x - a.x;
    double y2 = c.y - a.y;
    return x1 * y2 - x2 * y1;
}

bool cmp(Point a, Point b) {
    double c = cross(points[0], a, b);
    if (fabs(c) > EPS) {
        return c > 0;
    } else {
        return a.x < b.x;
    }
}

void convexHull() {
    int p0 = 0;
    for (int i = 1; i < n; i++) {
        if (points[i].y < points[p0].y || (points[i].y == points[p0].y && points[i].x < points[p0].x)) {
            p0 = i;
        }
    }
    swap(points[0], points[p0]);
    sort(points+1, points+n, cmp);
    stack<Point> st;
    st.push(points[0]);
    st.push(points[1]);
    for (int i = 2; i < n; i++) {
        while (st.size() >= 2) {
            Point p2 = st.top();
            st.pop();
            Point p1 = st.top();
            if (cross(p1, p2, points[i]) > 0) {
                st.push(p2);
                break;
            }
        }
        st.push(points[i]);
    }
    cout << st.size() << endl;
    while (!st.empty()) {
        Point p = st.top();
        st.pop();
        cout << p.x << " " << p.y << endl;
    }
}

int main() {
    cin >> n;
    for (int i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y;
        points[i] = Point(x, y);
    }
    convexHull();
    return 0;
}

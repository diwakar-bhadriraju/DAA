#include <iostream>
#include <cstdio>
#include <algorithm>
#include <set>
#include <cmath>

using namespace std;

const double EPS = 1e-9;
const int MAXN = 100005;

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

class Segment {
public:
    int p, q;
    Segment() {}
    Segment(int p, int q) : p(p), q(q) {}
};

int n;
Point points[2*MAXN];
Segment segments[MAXN];
set<pair<double, int> > eventSet;

double cross(Point a, Point b, Point c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool intersect(Segment s1, Segment s2) {
    double c1 = cross(points[s1.p], points[s1.q], points[s2.p]);
    double c2 = cross(points[s1.p], points[s1.q], points[s2.q]);
    double c3 = cross(points[s2.p], points[s2.q], points[s1.p]);
    double c4 = cross(points[s2.p], points[s2.q], points[s1.q]);
    if (((c1 > 0 && c2 < 0) || (c1 < 0 && c2 > 0)) && ((c3 > 0 && c4 < 0) || (c3 < 0 && c4 > 0))) {
        return true;
    }
    if (c1 == 0 && ((points[s2.p].x - points[s1.p].x) * (points[s1.q].x - points[s1.p].x) >= 0) && ((points[s2.p].y - points[s1.p].y) * (points[s1.q].y - points[s1.p].y) >= 0)) {
        return true;
    }
    if (c2 == 0 && ((points[s2.q].x - points[s1.p].x) * (points[s1.q].x - points[s1.p].x) >= 0) && ((points[s2.q].y - points[s1.p].y) * (points[s1.q].y - points[s1.p].y) >= 0)) {
        return true;
    }
    if (c3 == 0 && ((points[s1.p].x - points[s2.p].x) * (points[s2.q].x - points[s2.p].x) >= 0) && ((points[s1.p].y - points[s2.p].y) * (points[s2.q].y - points[s2.p].y) >= 0)) {
        return true;
    }
    if (c4 == 0 && ((points[s1.q].x - points[s2.p].x) * (points[s2.q].x - points[s2.p].x) >= 0) && ((points[s1.q].y - points[s2.p].y) * (points[s2.q].y - points[s2.p].y) >= 0)) {
        return true;
    }
    return false;
}

void addEvent(int i, int j) {
    if (i > j) {
        swap(i, j);
    }
    if (points[i].x > points[j].x) {
        eventSet.insert(make_pair(points[j].x, j));
    } else {
        eventSet.insert(make_pair(points[i].x, i));
    }
}

int main() {
    cin >> n;
    for (int i = 0; i < n; i++) {
        double x1, y1, x2, y2;
        cin >> x1 >> y1 >> x2 >> y2;
        points[2*i] = Point(x1, y1);
        points[2*i+1] = Point(x2, y2);
        segments[i] = Segment(2*i, 2*i+1);
        addEvent(2*i, 2*i+1);
    }
    int count = 0;
    while (!eventSet.empty()) {
        pair<double, int> currEvent = *eventSet.begin();
        eventSet.erase(eventSet.begin());
        int i = currEvent.second;
        if (i % 2 == 0) {
            intj = i + 1;
            for (set<pair<double, int> >::iterator it = eventSet.lower_bound(make_pair(points[i].x, -1)); it != eventSet.end(); it++) {
                int k = it->second;
                if (k == i || k == j) {
                    continue;
                }
                if (k % 2 == 0) {
                    if (intersect(segments[i/2], segments[k/2])) {
                        count++;
                    }
                } else {
                    if (intersect(segments[i/2], segments[(k-1)/2])) {
                        count++;
                    }
                }
            }
        } else {
            int j = i - 1;
            for (set<pair<double, int> >::iterator it = eventSet.lower_bound(make_pair(points[j].x, -1)); it != eventSet.end() && it->first <= points[i].x; it++) {
                int k = it->second;
                if (k == i || k == j) {
                    continue;
                }
                if (k % 2 == 0) {
                    if (intersect(segments[i/2], segments[k/2])) {
                        count++;
                    }
                } else {
                    if (intersect(segments[i/2], segments[(k-1)/2])) {
                        count++;
                    }
                }
            }
        }
        eventSet.erase(currEvent);
    }
    cout << count << endl;
    return 0;
}

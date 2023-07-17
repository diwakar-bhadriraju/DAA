#include<iostream>
#include<climits>
#include<algorithm>

using namespace std;

int assemblyLineScheduling(int entry[], int exit[], int a[][2], int t[][2], int n) {
    int f1[n], f2[n];
    f1[0] = entry[0] + a[0][0];
    f2[0] = entry[1] + a[0][1];
    for(int i=1; i<n; i++) {
        f1[i] = min(f1[i-1] + a[i][0], f2[i-1] + t[i-1][1] + a[i][0]);
        f2[i] = min(f2[i-1] + a[i][1], f1[i-1] + t[i-1][0] + a[i][1]);
    }
    return min(f1[n-1] + exit[0], f2[n-1] + exit[1]);
}

int main() {
    int entry[] = {2, 4};
    int exit[] = {3, 2};
    int a[][2] = {{7, 9}, {3, 4}, {4, 8}, {6, 5}, {8, 4}};
    int t[][2] = {{2, 3}, {4, 2}, {3, 1}, {1, 2}};
    int n = sizeof(a)/sizeof(a[0]);
    int minTime = assemblyLineScheduling(entry, exit, a, t, n);
    cout << "Minimum time to exit the assembly line = " << minTime << endl;
    return 0;
}

#include<iostream>
#include<cstring>
#include<algorithm>

using namespace std;

void printItemsIncluded(int K[][100], int wt[], int n, int W) {
    int i = n, w = W;
    cout << "Items included: ";
    while(i > 0 && w > 0) {
        if(K[i][w] != K[i-1][w]) {
            cout << i << " ";
            w -= wt[i-1];
        }
        i--;
    }
    cout << endl;
}

int knapSack(int W, int wt[], int val[], int n) {
    int K[n+1][100];
    for(int i=0; i<=n; i++) {
        for(int w=0; w<=W; w++) {
            if(i == 0 || w == 0) {
                K[i][w] = 0;
            } else if(wt[i-1] <= w) {
                K[i][w] = max(val[i-1] + K[i-1][w-wt[i-1]], K[i-1][w]);
            } else {
                K[i][w] = K[i-1][w];
            }
        }
    }
    printItemsIncluded(K, wt, n, W);
    return K[n][W];
}

int main() {
    int val[] = {60, 100, 120};
    int wt[] = {10, 20, 30};
    int W = 50;
    int n = sizeof(val)/sizeof(val[0]);
    cout << "Maximum value that can be obtained = " << knapSack(W, wt, val, n) << endl;
    return 0;
}

#include<iostream>
#include<climits>

using namespace std;

void printOptimalParenthesis(int s[][100], int i, int j) {
    if(i == j) {
        cout << "A" << i;
    } else {
        cout << "(";
        printOptimalParenthesis(s, i, s[i][j]);
        printOptimalParenthesis(s, s[i][j] + 1, j);
        cout << ")";
    }
}

void matrixChainMultiplication(int p[], int n) {
    int m[n][n], s[n][n];
    for(int i=1; i<n; i++) {
        m[i][i] = 0;
    }
    for(int L=2; L<n; L++) {
        for(int i=1; i<=n-L+1; i++) {
            int j = i+L-1;
            m[i][j] = INT_MAX;
            for(int k=i; k<j; k++) {
                int q = m[i][k] + m[k+1][j] + p[i-1] * p[k] * p[j];
                if(q < m[i][j]) {
                    m[i][j] = q;
                    s[i][j] = k;
                }
            }
        }
    }
    cout << "Minimum number of multiplications = " << m[1][n-1] << endl;
    cout << "Maximum number of multiplications = " << m[1][n-1] - p[0] * p[n-1] << endl;
    cout << "Optimal paranthesis: ";
    printOptimalParenthesis(s, 1, n-1);
}

int main() {
    int p[] = {30, 35, 15, 5, 10, 20, 25};
    int n = sizeof(p)/sizeof(p[0]);
    matrixChainMultiplication(p, n);
    return 0;
}

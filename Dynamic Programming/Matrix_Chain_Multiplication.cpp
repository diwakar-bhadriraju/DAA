#include <iostream>
#include <climits>
using namespace std;

void matrixChainOrder(int p[], int n, int m[][100], int s[][100], int maxM[][100], int minM[][100]) {
    for(int i=1; i<=n; i++) {
        m[i][i] = 0;
        maxM[i][i] = 0;
        minM[i][i] = 0;
    }
    for(int l=2; l<=n; l++) {
        for(int i=1; i<=n-l+1; i++) {
            int j = i + l - 1;
            m[i][j] = INT_MAX;
            maxM[i][j] = INT_MIN;
            minM[i][j] = INT_MAX;
            for(int k=i; k<=j-1; k++) {
                int q = m[i][k] + m[k+1][j] + p[i-1]*p[k]*p[j];
                if(q < m[i][j]) {
                    m[i][j] = q;
                    s[i][j] = k;
                }
                if(q > maxM[i][j]) {
                    maxM[i][j] = q;
                }
                if(q < minM[i][j]) {
                    minM[i][j] = q;
                }
            }
        }
    }
}

void printOptimalParenthesis(int s[][100], int i, int j, int n) {
    if(i == j) {
        cout << "A" << i;
    } else {
        cout << "(";
        printOptimalParenthesis(s, i, s[i][j], n);
        printOptimalParenthesis(s, s[i][j]+1, j, n);
        cout << ")";
    }
}

int main() {
    int n;
    cout << "Enter the number of matrices: ";
    cin >> n;
    int p[n+1];
    cout << "Enter the dimensions of the matrices: ";
    for(int i=0; i<=n; i++) {
        cin >> p[i];
    }
    int m[n+1][100], s[n+1][100], maxM[n+1][100], minM[n+1][100];
    matrixChainOrder(p, n, m, s, maxM, minM);
    cout << "Minimum number of scalar multiplications: " << m[1][n] << endl;
    cout << "Maximum number of scalar multiplications: " << maxM[1][n] << endl;
    cout << "Optimal Parenthesization: ";
    printOptimalParenthesis(s, 1, n, n);
    cout << endl;
    return 0;
}

#include<iostream>
#include<cstring>

using namespace std;

void printLCS(char* X, char* Y, int m, int n, int** L) {
    char LCS[100];
    int index = L[m][n];
    LCS[index] = '\0';
    int i = m, j = n;
    while(i > 0 && j > 0) {
        if(X[i-1] == Y[j-1]) {
            LCS[index-1] = X[i-1];
            i--;
            j--;
            index--;
        } else if(L[i-1][j] > L[i][j-1]) {
            i--;
        } else {
            j--;
        }
    }
    cout << "Length of LCS = " << L[m][n] << endl;
    cout << "LCS string = " << LCS << endl;
    for(int i=0; i<=m; i++) {
        delete[] L[i];
    }
    delete[] L;
}

void longestCommonSubsequence(char* X, char* Y, int m, int n) {
    int** L = new int*[m+1];
    for(int i=0; i<=m; i++) {
        L[i] = new int[n+1];
    }
    for(int i=0; i<=m; i++) {
        for(int j=0; j<=n; j++) {
            if(i == 0 || j == 0) {
                L[i][j] = 0;
            } else if(X[i-1] == Y[j-1]) {
                L[i][j] = L[i-1][j-1] + 1;
            } else {
                L[i][j] = max(L[i-1][j], L[i][j-1]);
            }
        }
    }
    printLCS(X, Y, m, n, L);
}

int main() {
    char X[] = "AGGTAB";
    char Y[] = "GXTXAYB";
    int m = strlen(X);
    int n = strlen(Y);
    longestCommonSubsequence(X, Y, m, n);
    return 0;
}

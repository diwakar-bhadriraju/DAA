#include<iostream>
#include<cmath>
#include <cstring>
using namespace std;

void printBoard(int** board, int n) {
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            cout << board[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

bool isSafe(int** board, int row, int col, int n) {
    for(int i=0; i<row; i++) {
        if(board[i][col] == 1) {
            return false;
        }
    }
    for(int i=row, j=col; i>=0 && j>=0; i--, j--) {
        if(board[i][j] == 1) {
            return false;
        }
    }
    for(int i=row, j=col; i>=0 && j<n; i--, j++) {
        if(board[i][j] == 1) {
            return false;
        }
    }
    return true;
}

void solveNQueens(int** board, int row, int n, int& count) {
    if(row == n) {
        count++;
        printBoard(board, n);
        return;
    }
    for(int col=0; col<n; col++) {
        if(isSafe(board, row, col, n)) {
            board[row][col] = 1;
            solveNQueens(board, row+1, n, count);
            board[row][col] = 0;
        }
    }
}

int countNQueens(int n) {
    int** board = new int*[n];
    for(int i=0; i<n; i++) {
        board[i] = new int[n];
        memset(board[i], 0, n*sizeof(int));
    }
    int count = 0;
    solveNQueens(board, 0, n, count);
    for(int i=0; i<n; i++) {
        delete[] board[i];
    }
    delete[] board;
    return count;
}

int main() {
    int n = 8;
    int count = countNQueens(n);
    cout << "Number of possible solutions = " << count << endl;
    return 0;
}

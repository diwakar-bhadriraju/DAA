#include<iostream>
#include<cstring>

using namespace std;

bool isSubsetSumBacktracking(int* set, int n, int sum) {
    if(sum == 0) {
        return true;
    }
    if(n == 0 && sum != 0) {
        return false;
    }
    if(set[n-1] > sum) {
        return isSubsetSumBacktracking(set, n-1, sum);
    }
    return isSubsetSumBacktracking(set, n-1, sum) || isSubsetSumBacktracking(set, n-1, sum-set[n-1]);
}

bool isSubsetSumRecursion(int* set, int n, int sum) {
    if(sum == 0) {
        return true;
    }
    if(n == 0 && sum != 0) {
        return false;
    }
    if(set[n-1] > sum) {
        return isSubsetSumRecursion(set, n-1, sum);
    }
    return isSubsetSumRecursion(set, n-1, sum) || isSubsetSumRecursion(set, n-1, sum-set[n-1]);
}

bool isSubsetSumDP(int* set, int n, int sum) {
    bool dp[n+1][sum+1];
    memset(dp, 0, sizeof(dp));
    for(int i=0; i<=n; i++) {
        dp[i][0] = true;
    }
    for(int i=1; i<=n; i++) {
        for(int j=1; j<=sum; j++) {
            if(set[i-1] > j) {
                dp[i][j] = dp[i-1][j];
            } else {
                dp[i][j] = dp[i-1][j] || dp[i-1][j-set[i-1]];
            }
        }
    }
    return dp[n][sum];
}

int main() {
    int set[] = {3, 34, 4, 12, 5, 2};
    int sum = 9;
    int n = sizeof(set)/sizeof(set[0]);
    cout << "Subset sum exists (backtracking): " << isSubsetSumBacktracking(set, n, sum) << endl;
    cout << "Subset sum exists (recursion): " << isSubsetSumRecursion(set, n, sum) << endl;
    cout << "Subset sum exists (dynamic programming): " << isSubsetSumDP(set, n, sum) << endl;
    return 0;
}

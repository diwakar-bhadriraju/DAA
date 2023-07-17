#include<iostream>
#include<cstring>
#include<algorithm>

using namespace std;

struct Job {
    int start, finish, profit;
};

bool jobComparison(Job a, Job b) {
    return (a.finish < b.finish);
}

int jobSelectionBranchAndBound(Job* jobs, int n) {
    sort(jobs, jobs+n, jobComparison);
    int* dp = new int[n];
    memset(dp, 0, sizeof(dp));
    dp[0] = jobs[0].profit;
    for(int i=1; i<n; i++) {
        int includeProfit = jobs[i].profit;
        int latestNonConflict = -1;
        for(int j=i-1; j>=0; j--) {
            if(jobs[j].finish <= jobs[i].start) {
                latestNonConflict = j;
                break;
            }
        }
        if(latestNonConflict != -1) {
            includeProfit += dp[latestNonConflict];
        }
        dp[i] = max(includeProfit, dp[i-1]);
    }
    int result = dp[n-1];
    delete[] dp;
    return result;
}

int jobSelectionDP(Job* jobs, int n) {
    sort(jobs, jobs+n, jobComparison);
    int* dp = new int[n];
    memset(dp, 0, sizeof(dp));
    dp[0] = jobs[0].profit;
    for(int i=1; i<n; i++) {
        int includeProfit = jobs[i].profit;
        for(int j=i-1; j>=0; j--) {
            if(jobs[j].finish <= jobs[i].start) {
                includeProfit += dp[j];
                break;
            }
        }
        dp[i] = max(includeProfit, dp[i-1]);
    }
    int result = dp[n-1];
    delete[] dp;
    return result;
}

int jobSelectionRecursion(Job* jobs, int n, int index, int* memo) {
    if(index == n) {
        return 0;
    }
    if(memo[index] != -1) {
        return memo[index];
    }
    int includeProfit = jobs[index].profit;
    int latestNonConflict = -1;
    for(int i=index-1; i>=0; i--) {
        if(jobs[i].finish <= jobs[index].start) {
            latestNonConflict = i;
            break;
        }
    }
    if(latestNonConflict != -1) {
        includeProfit += jobSelectionRecursion(jobs, n, latestNonConflict, memo);
    }
    int excludeProfit = jobSelectionRecursion(jobs, n, index+1, memo);
    int result = max(includeProfit, excludeProfit);
    memo[index] = result;
    return result;
}

int jobSelectionBacktracking(Job* jobs, int n, int index, int profit, int finishTime) {
    if(index == n) {
        return profit;
    }
    int includeProfit = 0;
    if(jobs[index].start >= finishTime) {
        includeProfit = jobSelectionBacktracking(jobs, n, index+1, profit+jobs[index].profit, jobs[index].finish);
    }
    int excludeProfit = jobSelectionBacktracking(jobs, n, index+1, profit, finishTime);
    return max(includeProfit, excludeProfit);
}

int main() {
    Job jobs[] = {{1, 2, 50}, {3, 5, 20}, {6, 19, 100}, {2, 100, 200}};
    int n = sizeof(jobs)/sizeof(jobs[0]);
    cout << "Job selection profit (branch and bound): " << jobSelectionBranchAndBound(jobs, n) << endl;
    cout << "Job selection profit (dynamic programming): " << jobSelectionDP(jobs, n) << endl;
    int* memo = new int[n];
    memset(memo, -1, sizeof(memo));
    cout << "Job selection profit (recursion): " << jobSelectionRecursion(jobs, n, 0, memo) << endl;
    delete[] memo;
    cout << "Job selection profit (backtracking): " << jobSelectionBacktracking(jobs, n, 0, 0, 0) << endl;
    return 0;
}

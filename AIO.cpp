#### N QUEENS
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

### Subset Sum
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


### Job Selection
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


### Karatsuba
#include <iostream>
#include <string>

using namespace std;

string add(string num1, string num2) {
    int n1 = num1.length(), n2 = num2.length();
    if (n1 < n2) {
        swap(num1, num2);
        swap(n1, n2);
    }
    string sum = "";
    int carry = 0, j = n2 - 1;
    for (int i = n1 - 1; i >= 0; i--) {
        int digit = num1[i] - '0';
        if (j >= 0) {
            digit += num2[j--] - '0';
        }
        digit += carry;
        carry = digit / 10;
        digit %= 10;
        sum = to_string(digit) + sum;
    }
    if (carry > 0) {
        sum = to_string(carry) + sum;
    }
    return sum;
}

string subtract(string num1, string num2) {
    int n1 = num1.length(), n2 = num2.length();
    if (n1 < n2) {
        swap(num1, num2);
        swap(n1, n2);
    }
    string diff = "";
    int carry = 0, j = n2 - 1;
    for (int i = n1 - 1; i >= 0; i--) {
        int digit = num1[i] - '0';
        if (j >= 0) {
            digit -= num2[j--] - '0';
        }
        digit -= carry;
        if (digit < 0) {
            digit += 10;
            carry = 1;
        } else {
            carry = 0;
        }
        diff = to_string(digit) + diff;
    }
    while (diff.length() > 1 && diff[0] == '0') {
        diff.erase(0, 1);
    }
    return diff;
}

string multiply(string num1, string num2) {
    int n1 = num1.length(), n2 = num2.length();
    if (n1 < n2) {
        swap(num1, num2);
        swap(n1, n2);
    }
    if (n1 == 0 || n2 == 0) {
        return "0";
    }
    if (n1 == 1) {
        int digit1 = num1[0] - '0', carry = 0;
        string prod = "";
        for (int i = n2 - 1; i >= 0; i--) {
            int digit2 = num2[i] - '0';
            int product = digit1 * digit2 + carry;
            carry = product / 10;
            product %= 10;
            prod = to_string(product) + prod;
        }
        if (carry > 0) {
            prod = to_string(carry) + prod;
        }
        return prod;
    }
    int mid = n1 / 2;
    string num1L = num1.substr(0, mid);
    string num1R = num1.substr(mid);
    string num2L = (n2 <= mid) ? "0" : num2.substr(0, n2 - mid);
    string num2R = num2.substr(max(0, n2 - mid));
    string prod1 = multiply(num1L, num2L);
    string prod2 = multiply(add(num1L, num1R), add(num2L, num2R));
    string prod3 = multiply(num1R, num2R);
    string term1 = prod1 + string(2 * (n1 - mid), '0');
    string term2 = subtract(subtract(prod2, prod1), prod3) + string(n1 - mid, '0');
    return add(add(term1, term2), prod3);
}

int main() {
    string num1 = "3141592653589793238462643383279502884197169399375105820974944592";
    string num2 = "2718281828459045235360287471352662497757247093699959574966967627";
    string product = multiply(num1, num2);
    cout << "Product = " << product << endl;
    return 0;
}

### Max SUbarray Sum
#include <iostream>
#include <climits>

using namespace std;

int maxCrossingSubarray(int arr[], int low, int mid, int high) {
    int leftSum = INT_MIN, rightSum = INT_MIN, sum = 0;
    for (int i = mid; i >= low; i--) {
        sum += arr[i];
        leftSum = max(leftSum, sum);
    }
    sum = 0;
    for (int i = mid + 1; i <= high; i++) {
        sum += arr[i];
        rightSum = max(rightSum, sum);
    }
    return leftSum + rightSum;
}

int maxSubarraySum(int arr[], int low, int high) {
    if (low == high) {
        return arr[low];
    }
    int mid = (low + high) / 2;
    int leftSum = maxSubarraySum(arr, low, mid);
    int rightSum = maxSubarraySum(arr, mid + 1, high);
    int crossSum = maxCrossingSubarray(arr, low, mid, high);
    return max(max(leftSum, rightSum), crossSum);
}

int main() {
    int arr[] = {-2, -3, 4, -1, -2, 1, 5, -3};
    int n = sizeof(arr) / sizeof(arr[0]);
    int maxSum = maxSubarraySum(arr, 0, n - 1);
    cout << "Maximum subarray sum = " << maxSum << endl;
    return 0;
}


### 0-1 knapsack
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


### Assembly Line Scheduling
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


### LCS
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


### MCM
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


### Convex Hull
#include <iostream>
#include <stack>
#include <algorithm>

struct Point {
    int x, y;
    Point(int a, int b) : x(a), y(b) {}
};

bool compare(Point p1, Point p2) {
    return p1.x < p2.x;
}
int orientation(Point p, Point q, Point r) {
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);
    if (val == 0) {
        return 0;
    }
    return (val > 0) ? 1 : 2;
}

void printConvexHull(Point points[], int n) {
    std::stack<Point> S;
    S.push(points[0]);
    S.push(points[1]);
    for (int i = 2; i < n; i++) {
        while (S.size() > 1 && orientation(S.top(), S.top(), points[i]) != 2) {
            S.pop();
        }
        S.push(points[i]);
    }
    while (!S.empty()) {
        Point p = S.top();
        std::cout << "(" << p.x << ", " << p.y << ")" << std::endl;
        S.pop();
    }
}

int main() {
    const int n = 3;
    Point points[n] = { {1, 2}, {3, 4}, {5, 6} };

    std::sort(points, points + n, compare);

    for (int i = 0; i < n; i++) {
        std::cout << "(" << points[i].x << ", " << points[i].y << ")" << std::endl;
    }

    return 0;
}


### Intersection

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


### Fractional Knapsack
#include <iostream>
#include <algorithm>

using namespace std;

struct Item {
    int weight;
    int value;
    double ratio;
};

bool compare(Item a, Item b) {
    return a.ratio > b.ratio;
}

double fractionalKnapsack(int capacity, Item items[], int n, bool included[]) {
    sort(items, items+n, compare);
    double totalValue = 0;
    for (int i = 0; i < n; i++) {
        if (capacity == 0) {
            return totalValue;
        }
        int weight = min(items[i].weight, capacity);
        totalValue += weight * items[i].ratio;
        capacity -= weight;
        included[i] = (weight == items[i].weight);
    }
    return totalValue;
}

int main() {
    int capacity, n;
    cout << "Enter the capacity of the knapsack: ";
    cin >> capacity;
    cout << "Enter the number of items: ";
    cin >> n;
    Item items[n];
    bool included[n];
    for (int i = 0; i < n; i++) {
        cout << "Enter weight and value of item " << i+1 << ": ";
        cin >> items[i].weight >> items[i].value;
        items[i].ratio = (double)items[i].value / items[i].weight;
        included[i] = false;
    }

    double maxValue = fractionalKnapsack(capacity, items, n, included);
    cout << "Maximum value that can be obtained = " << maxValue << endl;
    cout << "Items included: ";
    for (int i = 0; i < n; i++) {
        if (included[i]) {
            cout << i+1 << " ";
        }
    }
    cout << endl;

    return 0;
}


### Huffman Coding

#include <iostream>
#include <queue>
#include <unordered_map>

using namespace std;

class Node {
public:
    char ch;
    int freq;
    Node *left, *right;

    Node(char ch, int freq, Node *left = nullptr, Node *right = nullptr) {
        this->ch = ch;
        this->freq = freq;
        this->left = left;
        this->right = right;
    }

    ~Node() {
        delete left;
        delete right;
    }
};

struct Compare {
    bool operator() (Node* a, Node* b) {
        return a->freq > b->freq;
    }
};

void printCodes(Node* root, string code, string codes[]) {
    if (!root) {
        return;
    }
    if (root->ch != '\0') {
        codes[root->ch] = code;
        cout << root->ch << ": " << code << endl;
    }
    printCodes(root->left, code + "0", codes);
    printCodes(root->right, code + "1", codes);
}

void huffmanCoding(string text) {
    int freq[256] = {0};
    for (char ch : text) {
        freq[ch]++;
    }

    priority_queue<Node*, vector<Node*>, Compare> pq;
    for (int i = 0; i < 256; i++) {
        if (freq[i] > 0) {
            Node* node = new Node(i, freq[i]);
            pq.push(node);
        }
    }

    while (pq.size() > 1) {
        Node* left = pq.top();
        pq.pop();
        Node* right = pq.top();
        pq.pop();
        Node* parent = new Node('\0', left->freq + right->freq, left, right);
        pq.push(parent);
    }

    string codes[256];
    Node* root = pq.top();
    printCodes(root, "", codes);

    delete root;
}

int main() {
    string text = "hello world";
    huffmanCoding(text);
    return 0;
}


### KMP

#include <iostream>
#include <string>
#include <cstring>

using namespace std;

const int MAXN = 100005;

int matches[MAXN];

void computeLPS(char* pattern, int* lps) {
    int m = strlen(pattern);
    int len = 0;
    int i = 1;
    lps[0] = 0;
    while (i < m) {
        if (pattern[i] == pattern[len]) {
            len++;
            lps[i] = len;
            i++;
        } else {
            if (len != 0) {
                len = lps[len-1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
}

void kmpStringMatching(char* text, char* pattern) {
    int n = strlen(text);
    int m = strlen(pattern);
    int* lps = new int[m];
    computeLPS(pattern, lps);
    int i = 0, j = 0;
    int count = 0;
    while (i < n) {
        if (pattern[j] == text[i]) {
            i++;
            j++;
        }
        if (j == m) {
            matches[count++] = i-j;
            j = lps[j-1];
        } else if (i < n && pattern[j] != text[i]) {
            if (j != 0) {
                j = lps[j-1];
            } else {
                i++;
            }
        }
    }
    delete[] lps;
    if (count == 0) {
        cout << "Pattern not found" << endl;
    } else {
        cout << "Pattern found " << count << " time(s) at starting locations: ";
        for (int i = 0; i < count; i++) {
            cout << matches[i];
            if (i < count-1) {
                cout << ", ";
            }
        }
        cout << endl;
    }
}

int main() {
    char text[MAXN], pattern[MAXN];
    cout << "Enter the text: ";
    cin.getline(text, MAXN);
    cout << "Enter the pattern: ";
    cin.getline(pattern, MAXN);
    kmpStringMatching(text, pattern);
    return 0;
}


### Rabin Karp

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

using namespace std;

const int MAXN = 100005;
const int PRIME = 101;

int matches[MAXN];

void rabinKarpStringMatching(string text, string pattern) {
    int n = text.size();
    int m = pattern.size();
    int hashText = 0, hashPattern = 0;
    int d = 1;
    for (int i = 0; i < m-1; i++) {
        d = (d * 256) % PRIME;
    }
    for (int i = 0; i < m; i++) {
        hashPattern = (hashPattern * 256 + pattern[i]) % PRIME;
        hashText = (hashText * 256 + text[i]) % PRIME;
    }
    int count = 0;
    for (int i = 0; i <= n - m; i++) {
        if (hashPattern == hashText) {
            int j;
            for (j = 0; j < m; j++) {
                if (text[i+j] != pattern[j]) {
                    break;
                }
            }
            if (j == m) {
                matches[count++] = i;
            }
        }
        if (i < n - m) {
            hashText = ((hashText - text[i] * d) * 256 + text[i+m]) % PRIME;
            if (hashText < 0) {
                hashText += PRIME;
            }
        }
    }
    if (count == 0) {
        cout << "Pattern not found" << endl;
    } else {
        cout << "Pattern found " << count << " time(s) at starting locations: ";
        for (int i = 0; i < count; i++) {
            cout << matches[i];
            if (i < count-1) {
                cout << ", ";
            }
        }
        cout << endl;
    }
}

int main() {
    string text, pattern;
    cout << "Enter the text: ";
    getline(cin, text);
    cout << "Enter the pattern: ";
    getline(cin, pattern);
    rabinKarpStringMatching(text, pattern);
    return 0;
}


### Naive String matching

#include <iostream>
#include <string>
#include <cstring>

using namespace std;

const int MAXN = 100005;

int matches[MAXN];

void naiveStringMatching(char* text, char* pattern) {
    int n = strlen(text);
    int m = strlen(pattern);
    int count = 0;
    for (int i = 0; i <= n - m; i++) {
        int j;
        for (j = 0; j < m; j++) {
            if (text[i+j] != pattern[j]) {
                break;
            }
        }
        if (j == m) {
            matches[count++] = i;
        }
    }
    if (count == 0) {
        cout << "Pattern not found" << endl;
    } else {
        cout << "Pattern found " << count << " time(s) at starting locations: ";
        for (int i = 0; i < count; i++) {
            cout << matches[i];
            if (i < count-1) {
                cout << ", ";
            }
        }
        cout << endl;
    }
}

int main() {
    char text[MAXN], pattern[MAXN];
    cout << "Enter the text: ";
    cin.getline(text, MAXN);
    cout << "Enter the pattern: ";
    cin.getline(pattern, MAXN);
    naiveStringMatching(text, pattern);
    return 0;
}
********************************

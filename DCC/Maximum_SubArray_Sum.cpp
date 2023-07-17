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

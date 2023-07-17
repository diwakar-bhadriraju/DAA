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

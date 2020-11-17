#include <iostream>
#include "Solution.h"
using namespace std;
int main() {
    Solution solution;
    /*
    vector<Interval> arr = {{1,10}, {2, 9}, {3, 8}, {4, 7}, {5, 6}, {6, 6}};
    arr = solution.merge(arr);
    for(auto it:arr){
        cout<<it.start<<" "<<it.end<<endl;
    }*/
    int n = 2;
    vector<vector<int>> mat = solution.generateMatrix(n);
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            cout<<mat[i][j]<<" ";
        }
        cout<<endl;
    }
    return 0;
}

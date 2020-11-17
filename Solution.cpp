//
// Created by Babnish Vyas on 11/10/20.
//

#include "Solution.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <stack>
#include <queue>
#include <set>
#include <math.h>
using namespace std;

int Solution::maxSubArray(const vector<int> &A) {

    int current_max = A[0];
    int max_so_far = A[0];
    for(int i=1;i<A.size();++i){
        current_max = max(current_max+A[i],A[i]);
        max_so_far = max(current_max,max_so_far);
    }
    return max(max_so_far,current_max);
}

int Solution::coverPoints(vector<int> &X, vector<int> &Y) {
    int ans = 0;
    for(int i=1;i<X.size();++i){
        int x = X[i-1];
        int y = Y[i-1];
        int nx = X[i];
        int ny = Y[i];
        ans += max(abs(nx-x),abs(ny-y));
    }
    return ans;
}

int Solution::maxArr(vector<int> &A) {
    int max1 = INT_MIN; // (ai + i)
    int min1 = INT_MAX; // (ai + i)
    int max2 = INT_MIN; // (ai - i)
    int min2 = INT_MAX; // (ai - i)
    for(int i=0;i<A.size();++i){
        max1 = max(max1,A[i]+i);
        max2 = max(max2,A[i]-i);
        min1 = min(min1,A[i]+i);
        min2 = min(min2,A[i]-i);
    }
    return max(max1-min1,max2-min2);
}

int Solution::solve(int sum, vector<int> &arr) {
    unordered_map<long,int> mmap;
    long int _sum = 0;
    for(int i=0;i<arr.size();++i){
        _sum += arr[i];
    }
    if(_sum%3!=0){
        return 0;
    }
    long n3 = _sum/3;
    long n23 = 2*_sum/3;
    long current = 0;
    int ans = 0;
    for(int i=0;i<arr.size();++i){
        current += arr[i];
        if(current==n23){
            ans += mmap[n3];
        }
        mmap[current] += 1;
    }
    //for(auto it=mmap.begin();it!=mmap.end();++it){
    //    cout<<it->first<<" "<<it->second<<endl;
    //}
    if(_sum==0){
        int count = mmap[0];
        if(count<2){
            return 0;
        }
        int res = 0;
        for(int i=count-2;i>=1;i--){
            res += i;
        }
        return res;
    }
    return ans;
}

vector<int> Solution::flip(string s) {
    int max_sum = 0;
    int current_sum = 0;
    vector<int> ans(2,INT_MAX);
    int start = 0;
    for(int i=0;i<s.length();++i){
        int current;
        if(s[i]=='0'){
            current = 1;
        }else{
            current = -1;
        }
        current_sum += current;
        if(current_sum<0){
            start = i + 1;
            current_sum = 0;
        }
        if(current_sum>max_sum){
            max_sum = current_sum;
            ans[0] = start+1;
            ans[1] = i+1;
        }
    }
    if(ans[0]==INT_MAX) {
        return {};
    }
    return ans;
}

string Solution::largestNumber(const vector<int> &A) {
    vector<string> arr(A.size());
    for(int i=0;i<A.size();++i){
        arr[i] = to_string(A[i]);
    }
    sort(arr.begin(),arr.end(),[](const auto& a,const auto& b){
        string ab = a;
        ab.append(b);
        string ba = b;
        ba.append(a);
        return ab>ba;
    });
    string ans="";
    for(int i=0;i<arr.size();++i){
        ans.append(arr[i]);
    }
    int i = 0;
    while(i<ans.length()&&ans[i]=='0'){
        i++;
    }
    if(i==ans.length()){
        return "0";
    }else{
        return ans.substr(i,ans.length()-i+1);
    }
    return ans;
}


void Solution::rotate(vector<vector<int> > &mat) {
    int k = 0;
    int n = mat[0].size();
    //Transpose...........
    for(int i=0;i<=k;++i){
        if(k==n){
            break;
        }
        for(int j=0;j<=i;j++){
            swap(mat[i][j],mat[j][i]);
        }
        k++;
    }
    //Rotation around vertical axis............
    for(int i=0;i<mat.size();++i){
        for(int j=0;j<(n+1)/2;j++){
            swap(mat[i][j],mat[i][n-j-1]);
        }
    }
}

void Solution::setZeroes(vector<vector<int> > &A) {
    set<int> ii,jj;
    for(int i=0;i<A.size();++i){
        for(int j=0;j<A[0].size();++j){
            if(A[i][j]==0){
                ii.insert(i);
                jj.insert(j);
            }
        }
    }
    if(ii.size()==A.size()||jj.size()==A[0].size()){
        A.assign(A.size(),vector<int>(A[0].size(),0));
        return;
    }
    for(auto it=ii.begin();it!=ii.end();++it){
        int index = *it;
        for(int j=0;j<A[0].size();++j){
            A[index][j] = 0;
        }
    }

    for(auto it=jj.begin();it!=jj.end();++it){
        int index = *it;
        for(int i=0;i<A.size();++i){
            A[i][index] = 0;
        }
    }
}

vector<vector<int> > Solution::generateMatrix(int n) {
    vector<vector<int> > mat;
    mat.assign(n,vector<int>(n,0));
    int i_start = 0;
    int i_end = n-1;
    int direction = 0;
    int j_start = 0;
    int j_end = n-1;
    int index = 1;
    int i = i_start;
    int j = j_start;
    i_start++;
    while(index<=n*n){
        mat[i][j] = index;
        index++;
        if(direction==0){
            if(j==j_end){
                direction = 1;
                i++;
                j_end--;
            }else{
                j++;
            }
        }else if(direction==1){
            if(i==i_end){
                direction = 2;
                j--;
                i_end--;
            }else{
                i++;
            }
        }else if(direction==2){
            if(j==j_start){
                direction = 3;
                i--;
                j_start++;
            }else{
                j--;
            }
        }else if(direction==3){
            if(i==i_start){
                j = j_start;
                direction = 0;
                i_start++;
            }else{
                i--;
            }
        }

    }
    return mat;
}


vector<int> Solution::nextPermutation(vector<int> &arr) {
    int index = arr.size()-1;
    while(index>0&&arr[index-1]>=arr[index]){
        index--;
    }
    if(index==0){
        sort(arr.begin(),arr.end());
        return  arr;
    }
    int current = arr[index-1];
    int sindex = index;
    for(int j=index;j<arr.size();++j){
        if(arr[j]>current&&arr[j]<arr[sindex]){
            sindex = j;
        }
    }
    swap(arr[index-1],arr[sindex]);
    sort(arr.begin()+index,arr.end());
    return arr;
}

vector<int> Solution::findPerm(const string s, int n) {
    vector<int> arr(n);
    for(int i=1;i<=n;++i){
        arr[i-1] = i;
    }
    for(int i=1;i<=s.length();){
        if(s[i-1]=='I'){
            i++;
        }else{
            int start = i-1;
            while(i<s.length()&&s[i]=='D'){
                i++;
            }
            int end = i;
            i++;
            while(start<end){
                swap(arr[start],arr[end]);
                start++;
                end--;
            }
        }
    }
    return arr;
}

vector<Interval> Solution::insert(vector<Interval> &intervals, Interval newInterval){
    intervals.push_back(newInterval);
    sort(intervals.begin(),intervals.end(),[](const Interval& a,const Interval& b){
        return a.start<b.start;
    });
    Interval current = intervals[0];
    vector<Interval> res;
    for(int i=1;i<intervals.size();++i){
        if(current.end<=intervals[i].start){
            //not overlapping.......
            res.push_back(current);
            current = intervals[i];
            continue;
        }else{
            current.start = min(current.start,intervals[i].start);
            current.end = max(current.end,intervals[i].end);
        }
    }
    res.push_back(current);
    return res;
}

vector<Interval> Solution::merge(vector<Interval> &arr){
    vector<Interval> res;
    stack<Interval> stk;
    sort(arr.begin(),arr.end(),[](const Interval& a,const Interval& b){
        return a.start<b.start;
    });
    stk.push(arr[0]);
    for(int i=1;i<arr.size();++i){
        if(arr[i].start<=stk.top().end){
            Interval interval = stk.top();
            stk.pop();
            interval.start = min(interval.start, arr[i].start);
            interval.end = max(interval.end, arr[i].end);
            stk.push(interval);
        }else{
            res.push_back(stk.top());
            stk.pop();
            stk.push(arr[i]);
        }
    }
    while(!stk.empty()){
        res.push_back(stk.top());
        stk.pop();
    }
    return res;
}

int Solution::firstMissingPositive(vector<int> &arr) {
    int n = arr.size();
    for(int i=0;i<arr.size();++i){
        if(arr[i]==i+1){
            arr[i] = INT_MIN;
            continue;
        }
        if(arr[i]>n||arr[i]<=0){
            continue;
        }
        int index = arr[i];
        while(index<=n&&index>0){
            int next = arr[index-1];
            arr[index-1] = INT_MIN;
            index = next;
        }
    }
    for(int i=0;i<arr.size();++i){
        if(arr[i]!=INT_MIN){
            return i+1;
        }
    }
    return n+1;
}

vector<int> Solution::repeatedNumber(const vector<int> &arr){
    int _xor = arr[0];
    for(int i=1;i<arr.size();++i){
        _xor ^= arr[i];
    }
    for(int i=1;i<=arr.size();++i){
        _xor ^= i;
    }
    int mask = (_xor)&(~(_xor-1));
    int a=0,b=0;
    for(int i=1;i<=arr.size();++i){
        if(mask&i){
            a ^= i;
        }else{
            b ^= i;
        }
    }
    for(int i=0;i<arr.size();++i){
        if(mask&arr[i]){
            a ^= arr[i];
        }else{
            b ^= arr[i];
        }
    }
    int repeated,missing;
    for(int i=0;i<arr.size();++i){
        if(a==arr[i]){
            repeated = a;
            missing = b;
            break;
        }else if(b==arr[i]){
            repeated = b;
            missing = a;
            break;
        }
    }
    return {repeated,missing};
}

vector<int> Solution::wave(vector<int> &arr){
    sort(arr.begin(),arr.end());
    int n = (arr.size()%2==0)?arr.size():(arr.size()-1);
    for(int i=0;i<n-1;i+=2){
        swap(arr[i],arr[i+1]);
    }
    return  arr;
}

int Solution::solve(vector<vector<int> > &arr) {
    sort(arr.begin(),arr.end(),[](const vector<int>& a,const vector<int>& b){
        return a[0]<b[0];
    });
    priority_queue<int,vector<int>,priority_queue_compare> pq;
    pq.push(arr[0][1]);
    for(int i=1;i<arr.size();++i){
        if(arr[i][0]>=pq.top()){
            pq.pop();
        }
        pq.push(arr[i][1]);
    }
    return pq.size();
}

int Solution::candy(vector<int> &arr) {
    int n = arr.size();
    vector<int> left_to_right(n);
    vector<int> right_to_left(n);
    left_to_right[0] = 1;
    right_to_left[n-1] = 1;
    for(int i=1;i<n;++i){
        if(arr[i]>arr[i-1]){
            left_to_right[i] = left_to_right[i-1] + 1;
        }else{
            left_to_right[i] = 1;
        }

        if(arr[n-1-i]>arr[n-i]){
            right_to_left[n-1-i] = right_to_left[n-i] + 1;
        }else{
            right_to_left[n-1-i] = 1;
        }
    }

    int ans = 0;
    for(int i=0;i<n;++i){
        ans += max(left_to_right[i],right_to_left[i]);
    }
    return ans;
}

string Solution::find_prefix(string& s1,string& s2){
    int i = 0;
    while(i<s1.length()&&i<s2.length()){
        if(s1[i]!=s2[i]){
            return s1.substr(0,i);
        }
        i++;
    }
    if(i==s1.length()){
        return s1;
    }
    return s2;
}
string Solution::longestCommonPrefix(vector<string> &A) {
    if(A.size()==1){
        return A[0];
    }
    string prefix = find_prefix(A[0],A[1]);
    for(int i=2;i<A.size();++i){
        prefix = find_prefix(A[i],prefix);
    }
    return prefix;
}

int Solution::strStr(const string text, const string pattern) {
    if(pattern.length()==0 || text.length() == 0){

        return -1;
    }
    vector<int> longestSuffixPrefix(pattern.length(),0);
    int i = 0;
    int j = i+1;
    while(j<pattern.length()){
        if(pattern[i]==pattern[j]){
            longestSuffixPrefix[j] = i + 1;
            i++;
            j++;
        }else{
            if(i>0){
                i = longestSuffixPrefix[i-1];
            }else{
                j++;
            }
        }
    }

    int index = 0;
    i = 0;
    while(index<text.length()){
        if(pattern[i]==text[index]){
            i++;
            index++;
        }
        if(i==pattern.length()){
            return index-i;
        }else{
            if(index<text.length()&&pattern[i]!=text[index]){
                if(i>0){
                    i = longestSuffixPrefix[i-1];
                }else{
                    index++;
                }
            }
        }
    }
    return -1;
}

TreeNode* Solution::invertTree(TreeNode* root) {
    if(root==NULL){
        return NULL;
    }
    invertTree(root->left);
    invertTree(root->right);
    swap(root->left,root->right);
    return root;
}

int Solution::find_max_Depth(TreeNode *root, int &length) {
    if(root==NULL){
        return 0;
    }
    int left = find_max_Depth(root->left,length);
    int right = find_max_Depth(root->right,length);
    length = max(left+1,max(right+1,length));
    return max(left,right)+1;
}
int Solution::maxDepth(TreeNode* root){
    int ans = 0;
    find_max_Depth(root,ans);
    return ans;
}

void Solution::find_min_Depth(TreeNode *root, int &length, int current) {
    if(root==NULL){
        return;
    }
    if(root->left==NULL&&root->right==NULL){
        length = min(length,current+1);
        return;
    }
    find_min_Depth(root->left,length,current+1);
    find_min_Depth(root->right,length,current+1);
}
int Solution::minDepth(TreeNode* root){
    int length = INT_MAX;
    find_min_Depth(root,length,0);
    return length;
}

int Solution::isSameTree(TreeNode* root1, TreeNode* root2) {
    if(root1==NULL&&root2==NULL){
        return true;
    }
    if((root1&&!root2)||(!root1&&root2)){
        return false;
    }
    if(root1->val==root2->val){
        bool left = isSameTree(root1->left,root2->left);
        bool right = isSameTree(root1->right,root2->right);
        return left&&right;
    }
    return false;
}
bool Solution::checkSymmetry(TreeNode* r1,TreeNode* r2){
    if(r1==NULL&&r2==NULL){
        return true;
    }
    if((r1&&!r2)||(r2&&!r1)){
        return false;
    }
    if(r1->val==r2->val){
        bool left = checkSymmetry(r1->left,r2->right);
        bool right = checkSymmetry(r1->right,r2->left);
        return left&&right;
    }
    return false;
}
int Solution::isSymmetric(TreeNode* A) {
    return checkSymmetry(A,A);
}

void Solution::insert_tree(Node* root,int val){
    Node* current = root;
    for(int i=30;i>=0;i--){
        int bit = (val>>i)&1;
        if(bit==0){
            if(current->left==NULL){
                current->left = new Node();
            }
            current = current->left;
        }else{
            if(current->right==NULL){
                current->right = new Node();
            }
            current = current->right;
        }

    }
}
int Solution::find_xor_max(Node* root,int val){
    Node* current = root;
    int current_xor = 0;
    int M = pow(2,30);
    for(int i=30;i>=0;i--){
        int bit = (val>>i)&1;
        if(bit==0){
            if(current->right){
                current_xor += M;
                current = current->right;
            }else{
                current = current->left;
            }
        }else{
            if(current->left){
                current_xor += M;
                current = current->left;
            }else{
                current = current->right;
            }
        }
        M /= 2;
    }
    return current_xor;
}
int Solution::solve(vector<int> &arr, vector<int> &brr) {
    Node* root = NULL;
    root = new Node();
    int ans = 0;
    for(int i=0;i<arr.size();++i){
        insert_tree(root,arr[i]);
    }
    for(int i=0;i<brr.size();++i){
        ans = max(ans,find_xor_max(root,brr[i]));
    }
    return ans;
}

TreeNode* Solution::solve(TreeNode* root) {
    /*PostOrder Traversal*/
    if(root==NULL){
        return NULL;
    }
    root->left = solve(root->left);
    root->right = solve(root->right);
    if(root->left==NULL&&root->right==NULL){
        return root;
    }
    if(root->left==NULL){
        TreeNode* tmp = root->right;
        delete root;
        return tmp;
    }
    if(root->right==NULL){
        TreeNode* tmp = root->left;
        delete root;
        return tmp;
    }
    return root;
}

void Solution::find_path(TreeNode* root,int val,vector<int>& path,vector<int>& ans){
    if(root==NULL){
        return;
    }
    if(root->val==val){
        path.push_back(root->val);
        for(int i=0;i<path.size();++i){
            ans.push_back(path[i]);
        }
        return;
    }
    path.push_back(root->val);
    find_path(root->left,val,path,ans);
    find_path(root->right,val,path,ans);
    path.pop_back();
}
vector<int> Solution::solve(TreeNode* root, int val) {
    vector<int> path;
    vector<int> ans;
    find_path(root,val,path,ans);
    return ans;
}

int Solution::Height_post(TreeNode* root,bool& ans){
    if(root==NULL){
        return 0;
    }
    int _left = Height_post(root->left,ans);
    int _right = Height_post(root->right,ans);
    if(abs(_left-_right)>1){
        ans = false;
    }
    return max(_left,_right)+1;
}
int Solution::isBalanced(TreeNode* root) {
    bool ans = true;
    Height_post(root,ans);
    return  ans;
}
int Solution::solve(vector<int> &preorder) {
    stack<int> stk;
    int n = preorder.size();
    int root = INT_MIN;
    for(int i=0;i<n;++i){
        if(preorder[i]<root){
            return 0;
        }
        while(!stk.empty()&&stk.top()<preorder[i]){
            root = stk.top();
            stk.pop();
        }
        stk.push(preorder[i]);
    }
    return 1;
}

void Solution::find_element(TreeNode* root,int& k,int B,int& ans){
    if(root==NULL){
        return;
    }
    find_element(root->left,k,B,ans);
    k++;
    if(k==B){
        ans = root->val;
    }
    find_element(root->right,k,B,ans);
}
int Solution::kthsmallest(TreeNode* A, int B) {
    int k = 0;
    int ans = INT_MAX;
    find_element(A,k,B,ans);
    return ans;
}


int Solution::t2Sum(TreeNode* root, int val) {
    stack<TreeNode*> stk1,stk2;
    TreeNode* current1 = root;
    TreeNode* current2 = root;
    while(true){
        while(current1){
            stk1.push(current1);
            current1 = current1->left;
        }
        while(current2){
            stk2.push(current2);
            current2 = current2->right;
        }
        current1 = stk1.top();
        current2 = stk2.top();
        if(current1==current2&&current1!=NULL){
            return 0;
        }else if(current1->val+current2->val==val){
            return 1;
        }else if(current1->val+current2->val>val){
            current1 = NULL;
            current2 = current2->left;
            stk2.pop();
        }else if(current1->val+current2->val<val){
            current2 = NULL;
            current1 = current1->right;
            stk1.pop();
        }
    }
    return 0;
}

vector<int> Solution::find_Cousins(TreeNode *root, int val) {
    vector<pair<TreeNode*,TreeNode*> > q,r;
    if(root->val==val){
        return {};
    }
    q.push_back({root,NULL});
    bool flag = false;
    TreeNode* his_parent = NULL;
    while(q.size()){
        int n = q.size();
        for(int i=0;i<n;++i){
            auto current = q[i];
            TreeNode* node = current.first;
            if(node->val==val){
                flag = true;
                his_parent = current.second;
                break;
            }else{
                if(node->left){
                    r.push_back({node->left,node});
                }
                if(node->right){
                    r.push_back({node->right,node});
                }
            }
        }
        if(flag){
            vector<int> ans;
            for(auto it:q){
                if(it.second!=his_parent){
                    ans.push_back(it.first->val);
                }
            }
            return ans;
        }else{
            swap(r,q);
            r.clear();
        }
    }
    return {};
}

vector<int> Solution::right_view(TreeNode* root){
    queue<TreeNode*> q;
    q.push(root);
    vector<int> ans;
    int _last = root->val;
    while(!q.empty()){
        int n = q.size();
        for(int i=0;i<n;++i){
            auto it = q.front();
            q.pop();
            if(it->left){
                q.push(it->left);
            }
            if(it->right){
                q.push(it->right);
            }
            _last = it->val;
        }
        ans.push_back(_last);
    }
    return ans;
}

vector<vector<int> > Solution::zigzagLevelOrder(TreeNode* root) {
    vector<vector<int> > ans;
    queue<TreeNode*> q;
    q.push(root);
    int level = 0;
    while(!q.empty()){
        int n = q.size();
        for(int i=0;i<n;++i){
            auto it = q.front();
            q.pop();
            if(it->left){
                q.push(it->left);
            }
            if(it->right){
                q.push(it->right);
            }
            if(ans.size()==level){
                ans.push_back(vector<int>());
            }
            ans[level].push_back(it->val);
        }
        if(level%2==1){
            reverse(ans[level].begin(),ans[level].end());
        }
        level++;
    }
    return ans;
}

void Solution::connect(TreeLinkNode* root) {
    queue<TreeLinkNode*> q;
    q.push(root);
    while(!q.empty()){
        int n = q.size();
        for(int i=0;i<n;++i){
            TreeLinkNode* it = q.front();
            q.pop();
            if(it->left){
                q.push(it->left);
            }
            if(it->right){
                q.push(it->right);
            }
            if(i!=n-1){
                TreeLinkNode *next_right = q.front();
                it->next = next_right;
            }else{
                it->next = NULL;
            }
        }
    }
}

vector<vector<int> > Solution::verticalOrderTraversal(TreeNode* root) {
    if(root==NULL){
        return {};
    }
    map<int,vector<int> > vertical;
    queue<pair<TreeNode*,int> > q;
    q.push({root,0});
    while(!q.empty()){
        auto it = q.front();
        q.pop();
        if(it.first->left){
            q.push({it.first->left,it.second-1});
        }
        if(it.first->right){
            q.push({it.first->right,it.second+1});
        }
        vertical[it.second].push_back(it.first->val);
    }
    vector<vector<int> > ans;
    for(auto it=vertical.begin();it!=vertical.end();++it){
        ans.push_back(it->second);
    }
    return ans;
}

vector<int> Solution::inorderTraversal(TreeNode* A) {
    stack<TreeNode*> stk;
    vector<int> ans;
    TreeNode* root = A;
    while(root||!stk.empty()){
        while(root){
            stk.push(root);
            root = root->left;
        }
        root = stk.top();
        stk.pop();
        ans.push_back(root->val);
        root = root->right;
    }
    return ans;
}

vector<int> Solution::preorderTraversal(TreeNode* root) {
    stack<TreeNode*> stk;
    stk.push(root);
    vector<int> ans;
    while(!stk.empty()){
        TreeNode* current = stk.top();
        stk.pop();
        if(current->right){
            stk.push(current->right);
        }
        if(current->left){
            stk.push(current->left);
        }
        ans.push_back(current->val);
    }
    return ans;
}

int Solution::findMedian(vector<vector<int> > &mat) {
    int _min = INT_MAX;
    int _max = INT_MIN;
    int n = mat.size();
    int m = mat[0].size();
    for(int i=0;i<n;++i){
        _min = min(_min,mat[i][0]);
        _max = max(_max,mat[i][m-1]);
    }
    int median = (n*m+1)/2;
    while(_min<_max){
        int mid = _min + (_max-_min)/2;
        int places = 0;
        for(int i=0;i<n;++i){
            places += upper_bound(mat[i].begin(),mat[i].end(),mid)-mat[i].begin();
        }
        if(places<median){
            _min = mid+1;
        }else{
            _max = mid;
        }
    }
    return _min;
}


int Solution::sqrt(int num) {
    if(num==0){
        return 0;
    }
    long int low = 1;
    long int high = num;
    long int ans = low;
    while(low<high){
        long int mid = low + (high-low)/2;
        if(mid*mid==num){
            ans = mid;
            break;
        }
        if(mid*mid>num){
            high = mid;
        }else if(mid*mid<num){
            ans = max(ans,mid);
            low = mid+1;
        }
    }
    if(ans==INT_MIN){
        return 0;
    }
    return ans;
}

bool binary_search_row(vector<vector<int> >& matrix,int target,int row,int low,int high){
    if(low<=high){
        int mid = low+(high-low)/2;
        if(matrix[row][mid]==target){
            return true;
        }else if(matrix[row][mid]<target){
            return binary_search_row(matrix,target,row,mid+1,high);
        }else{
            return binary_search_row(matrix,target,row,low,mid-1);
        }
    }
    return false;
}
bool binary_search_col(vector<vector<int> >& matrix,int target,int low,int high){
    int n = matrix[0].size();
    if(low<=high){
        int mid = low+(high-low)/2;
        if(matrix[mid][n-1]==target){
            return true;
        }else if(matrix[mid][n-1]>target){
            if(matrix[mid][0]==target){
                return true;
            }else if(matrix[mid][0]<target){
                return binary_search_row(matrix,target,mid,0,n-1);
            }else{
                return binary_search_col(matrix,target,low,mid-1);
            }
        }else if(matrix[mid][n-1]<target){
            return binary_search_col(matrix,target,mid+1,high);
        }
    }
    return false;
}
bool searchMatrix_(vector<vector<int>>& matrix, int target) {
    int m = matrix.size();
    if(m==0){
        return false;
    }
    int n = matrix[0].size();
    if(n==0){
        return false;
    }
    bool result = binary_search_col(matrix,target,0,m-1);
    return result;
}
int Solution::searchMatrix(vector<vector<int> > &mat, int x) {
    return searchMatrix_(mat,x);
}

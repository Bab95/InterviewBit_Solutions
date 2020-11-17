//
// Created by Babnish Vyas on 11/10/20.
//

#ifndef INTERVIEWBIT_SOLUTION_H
#define INTERVIEWBIT_SOLUTION_H
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

struct Interval {
    int start;
    int end;
    Interval() : start(0), end(0) {}
    Interval(int s, int e) : start(s), end(e) {}
};

struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;
    TreeNode(int x) : val(x), left(NULL), right(NULL) {}
};

struct TreeLinkNode {
    int val;
    TreeLinkNode *left, *right, *next;
    TreeLinkNode(int x) : val(x), left(NULL), right(NULL), next(NULL) {}
};

struct priority_queue_compare{
    //min heap.........
    bool operator()(const int& a,const int& b){
        return a>b;
    }
};

struct Node{
    Node* left;
    Node* right;
    Node():left(NULL),right(NULL){}
};
class Solution {
public:

    /*MIN STEPS IN INFINITE GRID*/
    int coverPoints(vector<int> &X, vector<int> &Y);

    /*MAX SUM CONTIGOUS SUBARRAY*/
    int maxSubArray(const vector<int>& A);

    /*MAXIMUM ABSOLUTE DIFFERENCE*/
    int maxArr(vector<int>& A);

    /**/
    int solve(int sum,vector<int>& arr);

    /*FLIP*/
    vector<int> flip(string A);

    /*LARGEST NUMBER*/
    string largestNumber(const vector<int> &A);

    /*ROTATE MATRIX */
    void rotate(vector<vector<int> > &mat);

    /*SET MATRIX ZERO*/
    void setZeroes(vector<vector<int> > &A);

    /*SPIRAL ORDER MATRIX*/
    vector<vector<int> > generateMatrix(int A);

    /*NEXT PERMUTATION*/
    vector<int> nextPermutation(vector<int> &A);

    /*FIND PERMUTATION*/
    vector<int> findPerm(const string A, int B);

    /*MERGE INTERVALS*/
    vector<Interval> insert(vector<Interval> &intervals, Interval newInterval);

    /*MERGE OVERLAPPING INTERVALS*/
    vector<Interval> merge(vector<Interval> &A);

    /*FIRST MISSING POSITIVE*/
    int firstMissingPositive(vector<int> &A);

    /*REPEATED AND MISSING NUMBERS*/
    vector<int> repeatedNumber(const vector<int> &A);

    /*WAVE*/
    vector<int> wave(vector<int> &A);

    /*MEETING ROOMS*/
    int solve(vector<vector<int> > &A);

    /*DISTRIBUTE CANDY*/
    int candy(vector<int> &A);

    /*LONGEST COMMON PREFIX*/
    string longestCommonPrefix(vector<string> &A);
    string find_prefix(string& s1,string& s2);

    /*IMPLEMENT STRSTR*/
    int strStr(const string A, const string text);
/*----------------------------------------------------------TREE------------------------------------------------------*/
    /*INVERT BINARY TREE*/
    TreeNode* invertTree(TreeNode* A);

    /*MAX DEPTH OF BINARY TREE*/
    int maxDepth(TreeNode* A);
    int find_max_Depth(TreeNode* root,int& length);

    /*MIN DEPTH OF BINARY TREE*/
    int minDepth(TreeNode* root);
    void find_min_Depth(TreeNode* root,int& length,int current);

    /*SAME TREE's*/
    int isSameTree(TreeNode* root1, TreeNode* root2);

    /*SYMMETRICAL TREES*/
    int isSymmetric(TreeNode* root);
    bool checkSymmetry(TreeNode* r1,TreeNode* r2);

    /*MAX XOR BETWEEN TWO ARRAYS*/
    int solve(vector<int> &A, vector<int> &B);
    void insert_tree(Node* root,int val);
    int find_xor_max(Node* root,int val);

    /*REMOVE HALF NODES*/
    TreeNode* solve(TreeNode* A);

    /*PATH TO A GIVEN NODE*/
    vector<int> solve(TreeNode* A, int B);
    void find_path(TreeNode* root,int val,vector<int>& path,vector<int>& ans);

    /*BALANCED TREE*/
    int isBalanced(TreeNode* A);
    int Height_post(TreeNode* root,bool& ans);

    /*VALID BST*/
    int solve(vector<int> &A);

    /*KTH SMALLEST ELEMENT*/
    int kthsmallest(TreeNode* A, int B);
    void find_element(TreeNode* root,int& k,int B,int& ans);

    /*2SUM*/
    int t2Sum(TreeNode* A, int B);

    /*BST ITERATOR*/


    /*COUSINS IN BINARY TREE*/
    vector<int> find_Cousins(TreeNode* A, int B);

    /*RIGHT VIEW OF BINARY TREE*/
    vector<int> right_view(TreeNode* A);

    /*ZIGZAG TRAVERSAL*/
    vector<vector<int> > zigzagLevelOrder(TreeNode* A);

    /*CONNECT NODES AT SAME LEVEL*/
    void connect(TreeLinkNode* A);

    /*VERTICALORDER TRAVERSAL*/
    vector<vector<int> > verticalOrderTraversal(TreeNode* A);

    /*INORDER TRAVERSAL*/
    vector<int> inorderTraversal(TreeNode* A);

    /*PREORDER TRAVERSAL*/
    vector<int> preorderTraversal(TreeNode* A);

    /*MEDIAN IN A MATRIX ROW WISE SORTED*/
    int findMedian(vector<vector<int> > &mat);

    /*SQRT OF A NUMBER*/
    int sqrt(int A);

    int searchMatrix(vector<vector<int> > &mat, int x);


};


#endif //INTERVIEWBIT_SOLUTION_H

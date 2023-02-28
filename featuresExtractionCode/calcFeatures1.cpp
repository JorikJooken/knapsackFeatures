// compile this file using: g++ -g -std=gnu++11 -O2 calcFeatures1.cpp -o calcFeatures1Executable -lquadmath

#include <bits/stdc++.h>
#include <quadmath.h>
//#include <boost/multiprecision/__float128.hpp>

using namespace std;
//using namespace boost::multiprecision;

vector<long long> profits;
vector<long long> weights;
vector<__float128> prefix_sums;
vector<__float128> prefix_sums_squares;

int n;
long long capacity;

__float128 inf=1e50Q;
__float128 eps=1e-9Q;

vector< vector<__float128> > dp;
vector< vector<int> > argmin;

// works in O(n)
void SMAWK(vector<int> &rows, vector<int> &columns, int dp_idx)
{
    if(rows.size()==1)
    {
        int row=rows[0];
        for(int col : columns)
        {
            __float128 cost;
            if(col>row) cost=inf; 
            else
            {
                __float128 avg=(prefix_sums[row]-prefix_sums[col-1])/(row-col+1);
                cost=dp[dp_idx-1][col-1]+
                     (prefix_sums_squares[row]-prefix_sums_squares[col-1])-
                     2*avg*(prefix_sums[row]-prefix_sums[col-1])+
                     (row-col+1)*avg*avg;
            }
            if(cost<0) cost=0; // for numerical reasons
            if(cost<dp[dp_idx][row]-eps)
            {
                dp[dp_idx][row]=cost;
                argmin[dp_idx][row]=col;
            }
        }
        return;   
    }
    if(rows.size()<columns.size())
    {
        vector<int> surviving_columns;
        int sz=0;
        for(int col : columns)
        {
            while(sz>0)
            {
                int row=rows[sz-1];
                __float128 cost1;
                if(surviving_columns.back()>row) cost1=inf; 
                else
                {
                    __float128 avg=(prefix_sums[row]-prefix_sums[surviving_columns.back()-1])/(row-surviving_columns.back()+1);
                    cost1=dp[dp_idx-1][surviving_columns.back()-1]+
                         (prefix_sums_squares[row]-prefix_sums_squares[surviving_columns.back()-1])-
                         2*avg*(prefix_sums[row]-prefix_sums[surviving_columns.back()-1])+
                         (row-surviving_columns.back()+1)*avg*avg;
                }
                __float128 cost2;
                if(col>row) cost2=inf; 
                else
                {
                    __float128 avg=(prefix_sums[row]-prefix_sums[col-1])/(row-col+1);
                    cost2=dp[dp_idx-1][col-1]+
                         (prefix_sums_squares[row]-prefix_sums_squares[col-1])-
                         2*avg*(prefix_sums[row]-prefix_sums[col-1])+
                         (row-col+1)*avg*avg;
                }
                if(cost1<0) cost1=0; // for numerical reasons
                if(cost2<0) cost2=0; // for numerical reasons
                if(cost1-eps>cost2)
                {
                    surviving_columns.pop_back();
                    sz--;
                }
                else break;
            }
            if(surviving_columns.size()<rows.size())
            {
                surviving_columns.push_back(col);
                sz++;
            }            
        }
        SMAWK(rows,surviving_columns,dp_idx);
        return;
    }
    else
    {
        vector<int> new_rows;
        for(int i=0; i<rows.size(); i+=2)
        {
            new_rows.push_back(rows[i]);
        }
        SMAWK(new_rows,columns,dp_idx);
        for(int i=1; i<rows.size(); i+=2)
        {
            int row=rows[i];
            int upper_bound=row;
            if(i+1<rows.size()) upper_bound=min(upper_bound,argmin[dp_idx][rows[i+1]]);
            for(int col=argmin[dp_idx][rows[i-1]]; col<=upper_bound; col++)
            {
                __float128 cost;
                if(col>row) cost=inf; 
                else
                {
                    __float128 avg=(prefix_sums[row]-prefix_sums[col-1])/(row-col+1);
                    cost=dp[dp_idx-1][col-1]+
                         (prefix_sums_squares[row]-prefix_sums_squares[col-1])-
                         2*avg*(prefix_sums[row]-prefix_sums[col-1])+
                         (row-col+1)*avg*avg;
                }
                if(cost<0) cost=0; // for numerical reasons

                if(cost<dp[dp_idx][row]-eps)
                {
                    dp[dp_idx][row]=cost;
                    argmin[dp_idx][row]=col;
                }
            }
        }
        return;
    }
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    cin >> n;
    // profits are unrelated to this problem in fact
    profits.assign(n,-1);
    weights.assign(n,-1);
    for(int i=0; i<n; i++)
    {
        long long id, pr, wei;
        cin >> id >> pr >> wei;
        profits[i]=pr;
        weights[i]=wei;
    }
    cin >> capacity;
    sort(weights.begin(),weights.end());
    prefix_sums.assign(n,0);
    prefix_sums_squares.assign(n,0);
    for(int i=0; i<n; i++)
    {
        if(i-1>=0)
        {
            prefix_sums[i]+=prefix_sums[i-1];
            prefix_sums_squares[i]+=prefix_sums_squares[i-1];
        }
        prefix_sums[i]+=weights[i];
        prefix_sums_squares[i]+=1.0Q*weights[i]*weights[i];
    }
    // dp[i][j] : min cost to make i+1 groups for items 0..j
    vector<__float128> dp_row(n,inf);
    dp.push_back(dp_row);
    vector<int> argmin_row(n,-1);
    argmin.push_back(argmin_row);
    // base case
    for(int i=0; i<n; i++)
    {
        __float128 avg=prefix_sums[i]/(i+1);
        // sum j=0..i (w_j-avg)^2 = (sum j=0..i w_j^2)-(2*avg*sum j=0..i w_j)+(i+1)*avg^2
        dp[0][i]=prefix_sums_squares[i]-2*avg*prefix_sums[i]+(i+1)*avg*avg;
        //dp[0][i]=max(dp[0][i],0.0); // for numerical reasons
    }    
    // recursive case
    for(int i=1; i<n; i++)
    {
        vector<__float128> new_dp_row(n,inf);
        dp.push_back(new_dp_row);
        vector<int> new_argmin_row(n,-1);
        argmin.push_back(new_argmin_row);
        
        // execute SMAWK algorithm
        // refers to the rows and columns of a virtual cost matrix that measures the cost
        vector<int> rows;
        for(int j=i;j<n;j++) rows.push_back(j);
        vector<int> columns;
        for(int j=i;j<n;j++) columns.push_back(j);
        // populate new row of dp table
        SMAWK(rows,columns,i);
        // later, put break condition here to find g*
        if(dp[i][n-1] > 0.9Q*dp[i-1][n-1]) break;
    }
    int g_star=dp.size()-1;
    int current_first_idx=g_star-1;
    int current_second_idx=n-1;
    vector<int> ends_of_groups_for_increasing_weights;
    while(current_first_idx>=0)
    {
        ends_of_groups_for_increasing_weights.push_back(current_second_idx);
        current_second_idx=argmin[current_first_idx][current_second_idx]-1;
        current_first_idx--;
    }
    vector<int> ends_of_groups_for_decreasing_weights;
    for(int i=0; i+1<ends_of_groups_for_increasing_weights.size(); i++)
    {
        int new_end=n-ends_of_groups_for_increasing_weights[i+1]-2;
        ends_of_groups_for_decreasing_weights.push_back(new_end);
    }
    ends_of_groups_for_decreasing_weights.push_back(n-1);

    printf("%d\n",g_star);
    char s[256];
    quadmath_snprintf(s, 256, "%Qf", dp[g_star-1][n-1]);
    printf("%s\n", s);
    for(int x : ends_of_groups_for_decreasing_weights)
    {
        printf("%d ", x);
    }
    printf("\n");
    /*
    for(int x : ends_of_groups_for_increasing_weights)
    {
        cerr << x << " ";
    }
    cerr << endl;*/

    
    /*for(int i=0; i<dp.size(); i++)
    {
        char s[256];
        quadmath_snprintf(s, 256, "%Qf", dp[i][n-1]);
        printf("%d -> %s\n", i, s);
    }
    //cout << dp.size() << endl;*/
    return 0;
}

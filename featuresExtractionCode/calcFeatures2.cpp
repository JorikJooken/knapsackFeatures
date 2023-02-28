// compile as g++ -g -std=c++11 -O2 calcFeatures2.cpp -o calcFeaturesExecutable
#include <bits/stdc++.h>

using namespace std;


vector<long long> profits;
vector<long long> weights;

int n;
long long capacity;

inline double added(double log2_a, double log2_b)
{
    if(log2_a<=-1000 && log2_b<=-1000) return -1000;
    if(log2_b>log2_a)
    {
        swap(log2_a,log2_b);
    }
    double diff=log2_a-log2_b;
    //a is much larger
    if(diff>1000.0)
    {
        return log2_a;
    }
    else
    {
        return log2_a+log2(1+pow(2.0,-diff));
    }
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    cin >> n;
    // profits are unrelated to whether a solution is inclusionwise maximal or not
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
    long long sum_first_i_minus_one_weights=0;
    for(int i=0;i<n;i++) sum_first_i_minus_one_weights+=weights[i];

    vector<double> log2_counts(capacity+1,-1000);
    log2_counts[0]=0.0; // pow(2.0,0.0)=1
    long long max_reached=0;
    double log2_numberInclusionwiseMaximalSolutions=-1000; //pow(2.0,-1000) is approximately equal to 0
    map<long long,double> all_log2_counts;
    // Time complexity: O(n*c)
    for(int i=n-1; i>=0; i--)
    {
        sum_first_i_minus_one_weights-=weights[i];
        // first i-1 items are selected, item i is not selected
        long long at_least=capacity+1-sum_first_i_minus_one_weights-weights[i];
        long long at_most=capacity-sum_first_i_minus_one_weights;
        if(i+1<n)
        {
            long long start=min(max_reached,capacity-weights[i+1]);
            for(long long j=start; j>=0; j--)
            {
                if(log2_counts[j]>-1000)
                {
                    log2_counts[j+weights[i+1]]=added(log2_counts[j+weights[i+1]],log2_counts[j]);
                    max_reached=max(max_reached,j+weights[i+1]);
                }
            }
        }
        long long start=max(0LL,at_least);
        long long end=min(capacity,at_most);
        for(long long j=start; j<=end; j++)
        {
            log2_numberInclusionwiseMaximalSolutions=added(log2_numberInclusionwiseMaximalSolutions,log2_counts[j]);
            double current_res=-1000;
            auto it=all_log2_counts.find(j+sum_first_i_minus_one_weights);
            if(it != all_log2_counts.end())
            {
                current_res=it->second;
            }
            if(current_res>-1000 || log2_counts[j]>-1000) all_log2_counts[j+sum_first_i_minus_one_weights]=added(current_res,log2_counts[j]);
        }
    }
    cout << "logarithm: " << fixed << setprecision(18) << log2_numberInclusionwiseMaximalSolutions << endl;
    cout << "pow(2,logarithm): " << fixed << setprecision(18) << pow(2.0,log2_numberInclusionwiseMaximalSolutions) << endl;
    cout << "Several statistics about the weight distribution of inclusionwise maximal solutions follow: " << endl;
    cout << "Number of different weights: " << all_log2_counts.size() << endl;
    pair<long long, double> smallest=make_pair(1e18,1e18);
    pair<long long, double> largest=make_pair(-1,-1);
    double log2_average=-1000.0;
    for(auto it=all_log2_counts.begin(); it != all_log2_counts.end(); it++)
    {
        smallest=min(smallest,make_pair(it->first,it->second));
        largest=max(largest,make_pair(it->first,it->second));
        log2_average=added(log2_average,log2((double)it->first)+it->second);
    }
    log2_average-=log2_numberInclusionwiseMaximalSolutions;
    cout << "Minimum weight and log(amount): " << smallest.first << " " << smallest.second << endl;
    cout << "Maximum weight and log(amount)= " << largest.first << " " << largest.second << endl;
    double average=pow(2.0,log2_average);
    cout << "Average weight: " << fixed << setprecision(18) << average << endl;
    double log2_sigma_squared=-1000.0;
    for(auto it=all_log2_counts.begin(); it != all_log2_counts.end(); it++)
    {
        if(fabs(it->first-average)<=1e-100) continue;
        log2_sigma_squared=added(log2_sigma_squared,it->second+log2(fabs(it->first-average))+log2(fabs(it->first-average)));
    }
    log2_sigma_squared-=log2_numberInclusionwiseMaximalSolutions;
    cout << "Sigma squared: " << fixed << setprecision(18) << pow(2.0, log2_sigma_squared) << endl;
    cout << "Percentiles: " << endl;
    double log2_cumulative_amount=-1000.0;
    auto it=all_log2_counts.begin();
    for(int i=1; i<=100; i++)
    {
        double target=log2(1.0*i/100)+log2_numberInclusionwiseMaximalSolutions;
        while(it != all_log2_counts.end() && log2_cumulative_amount<target)
        {
            log2_cumulative_amount=added(log2_cumulative_amount,it->second);
            it++;
        }
        it--;
        cout << i << "/100 percentile: " << it->first << endl;
        it++;        
    }
    /*
    cout << "All inclusionwise maximal solutions in format (weight_sum, log(amount)):" << endl;
    for(auto it=all_log2_counts.begin(); it != all_log2_counts.end(); it++)
    {
        cout << "(" << it->first << ", " << fixed << setprecision(18) << it->second << ")" << endl;
    }*/
    return 0;
}

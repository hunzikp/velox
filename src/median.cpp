#include <vector>
#include <algorithm>

using namespace std;

double median(vector<double> scores) {
  double median;
  size_t size = scores.size();

  sort(scores.begin(), scores.end());

  if (size  % 2 == 0) {
    median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
  }
  else {
    median = scores[size / 2];
  }

  return median;
}

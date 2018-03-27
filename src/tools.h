#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
    /**
     * Constructor.
     */
    Tools();

    /**
     * Destructor.
     */
    virtual ~Tools();

    /**
     * A helper method to calculate RMSE.
     */
    VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

    static double normalizeAngle(double angle) {
        double a = fmod(angle + M_PI, 2 * M_PI);
        return a >= 0 ? (a - M_PI) : (a + M_PI);
    }
};

#endif /* TOOLS_H_ */
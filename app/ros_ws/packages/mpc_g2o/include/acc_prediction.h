#ifndef ACC_PREDICTION_H
#define ACC_PREDICTION_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <deque>

class AccPrediction
{
public:
    std::vector<double> perform_v_p_prediction_task(int N, const double v_p, const double a_p, double a_p_weighted)
    {
        std::vector<double> F(N + 1);
        std::iota(F.begin(), F.end(), 0); // F = [0, 1, 2, ..., N]

        std::vector<double> v_p_prediction(F.size(), 0);
        if (a_p_weighted >= 0)
        {
            a_p_weighted = std::min(a_p, 0.5 * a_p_weighted);
            for (size_t i = 0; i < F.size(); ++i)
            {
                v_p_prediction[i] = a_p_weighted * F[i] + v_p;
                if (i >= 3)
                    v_p_prediction[i] = v_p_prediction[2]; // Fix values after the 3rd entry
            }
        }
        else
        {
            a_p_weighted = std::max({a_p, a_p_weighted, -0.2});
            for (size_t i = 0; i < F.size(); ++i)
            {
                v_p_prediction[i] = a_p_weighted * F[i] + v_p;
                if (v_p_prediction[i] < 0)
                    v_p_prediction[i] = 0;
            }
        }
        return v_p_prediction;
    }

    std::vector<double> perform_road_grad_prediction_task(int N, double gradient)
    {
        return std::vector<double>(N + 1, gradient);
    }

    std::vector<double> perform_v_min_prediction_task(int N, double v_min)
    {
        return std::vector<double>(N + 1, v_min);
    }

    std::vector<double> perform_v_max_prediction_task(int N, double v_max)
    {
        return std::vector<double>(N + 1, v_max);
    }

    std::vector<std::vector<double>> compute_mpc_prediction(int N, const double v_p, const double a_p,
                                                            double a_p_weighted,
                                                            double v_min, double v_max, double road_grad)
    {
        auto v_p_pred = perform_v_p_prediction_task(N, v_p, a_p, a_p_weighted);
        auto road_grad_pred = perform_road_grad_prediction_task(N, road_grad);
        auto v_min_pred = perform_v_min_prediction_task(N, v_min);
        auto v_max_pred = perform_v_max_prediction_task(N, v_max);

        return {v_p_pred, road_grad_pred, v_min_pred, v_max_pred};
    }
};

#endif // ACC_PREDIICTION_H

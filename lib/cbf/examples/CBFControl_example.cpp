//
// Created by lishuo on 9/2/24.
//

#include <model/DoubleIntegratorXYYaw.h>
#include <cbf/detail/cbf.h>
#include <cbf/controller/CBFControl.h>
#include <nlohmann/json.hpp>
#include <fstream>

constexpr unsigned int DIM = 3U;
using State = model::State<double, DIM>;
using VectorDIM = math::VectorDIM<double, DIM>;

math::VectorDIM<double, DIM> criticallyDampedSpringControl(const State& current_state, const VectorDIM& target, const double spring_constant) {
    VectorDIM Fs;
    Fs = spring_constant * target - current_state.pos_;
    VectorDIM Fd;
    Fd = -current_state.vel_ * 2 * sqrt(spring_constant);
    return Fs + Fd;
}

int main() {
    using FovCBF = cbf::FovCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using json = nlohmann::json;

    using Vector = math::Vector<double>;
    // load experiment config
    std::string experiment_config_filename = "../../../config/config.json";
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);
    double h = experiment_config_json["mpc_params"]["h"];

    double fov_beta = experiment_config_json["fov_cbf_params"]["beta"];
    double fov_Ds = experiment_config_json["fov_cbf_params"]["Ds"];
    double fov_Rs = experiment_config_json["fov_cbf_params"]["Rs"];


    // json for record
    std::string JSON_FILENAME = "../../../tools/CBFXYYawStates.json";
    json states;
    states["dt"] = h;
    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(h);
    // init cbf
    std::shared_ptr<FovCBF> fov_cbf = std::make_unique<FovCBF>(fov_beta, fov_Ds, fov_Rs);

    // load the tasks
    std::vector<State> init_states;
    std::vector<Vector> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    for (size_t i = 0; i < num_robots; ++i) {
        // load init states
        State init_state;
        init_state.pos_ << so_json[i][0], so_json[i][1], so_json[i][2];
        init_state.vel_ << VectorDIM::Zero();
        init_states.push_back(init_state);
        Vector current_state(6);
        current_state << init_state.pos_, init_state.vel_;
        current_states.push_back(current_state);
        // load target positions
        VectorDIM target_pos;
        target_pos << sf_json[i][0], sf_json[i][1], sf_json[i][2];
        target_positions.push_back(target_pos);
    }

    // control loop
    int loop_idx = 0;
    while (loop_idx < 200) {
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            VectorDIM the_other_robot_position;
            if (robot_idx == 0) {
                the_other_robot_position = init_states.at(1).pos_;
            } else {
                the_other_robot_position = init_states.at(0).pos_;
            }

            // compute the desired control
            const VectorDIM& target_pos = target_positions.at(robot_idx);
            const VectorDIM& current_pos = init_states.at(robot_idx).pos_;
            VectorDIM desired_u = criticallyDampedSpringControl(init_states.at(robot_idx), target_pos, 1.);


            State next_init_state = pred_model_ptr->applyInput(init_states.at(robot_idx), desired_u);
            init_states.at(robot_idx) = next_init_state;

            current_states.at(robot_idx)(0) = init_states.at(robot_idx).pos_(0);
            current_states.at(robot_idx)(1) = init_states.at(robot_idx).pos_(1);
            current_states.at(robot_idx)(2) = init_states.at(robot_idx).pos_(2);
            current_states.at(robot_idx)(3) = init_states.at(robot_idx).vel_(0);
            current_states.at(robot_idx)(4) = init_states.at(robot_idx).vel_(1);
            current_states.at(robot_idx)(5) = init_states.at(robot_idx).vel_(2);

            // log the state
            Vector x_t(6);
            x_t = current_states.at(robot_idx);
            states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
        }
        loop_idx += 1;
    }

    // write states to json
    std::cout << "writing to file " << JSON_FILENAME << ".\n";
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;

    return 0;
}
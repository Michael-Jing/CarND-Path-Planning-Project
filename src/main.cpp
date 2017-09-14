#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y) {
    int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];

    double next_map_x = maps_x[closestWaypoint + 1];
    double next_map_y = maps_y[closestWaypoint + 1];
    double next_cur_dif_x = next_map_x - map_x;
    double next_cur_dif_y = next_map_y - map_y;
    double ego_map_dif_x = x - map_x;
    double ego_map_dif_y = y - map_y;
    if ((next_cur_dif_x * ego_map_dif_x + next_cur_dif_y * ego_map_dif_y) >= 0) {
        return closestWaypoint + 1;
    } else {
        return closestWaypoint;
    }


}



vector<double> interopolate_normal_vector(vector<double> v1, vector<double> v2, double w1, double w2) {
    double x = v1[0] * w1 + v2[0] * w2;
    double y = v1[1] * w1 + v2[1] * w2;
    double length = sqrt(x * x + y * y);
    return {x / length, y / length};

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, vector<double> maps_x, vector<double> maps_y, vector<double> maps_s,
        vector<double> maps_dx,
    vector<double> maps_dy)
{
	int next_wp = NextWaypoint(x, y, maps_x, maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0) // the road is looped in this project.
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto T
    double wp_distance = sqrt(n_x * n_x + n_y * n_y);
    double s_to_prev_wp = (x_x * n_x + x_y * n_y) / wp_distance;
    double s_to_next_wp = wp_distance - s_to_prev_wp;

    double frenet_s = maps_s[prev_wp] + s_to_prev_wp;

    vector<double> normal_vector = interopolate_normal_vector({maps_dx[prev_wp], maps_dy[prev_wp]},
                                                              {maps_dx[next_wp], maps_dy[next_wp]},
                                                                 s_to_next_wp, s_to_prev_wp);

	double frenet_d = x_x * normal_vector[0] + x_y * normal_vector[1];

	return {frenet_s,frenet_d};

}


// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_x, vector<double> maps_y, vector<double> maps_s,
                     vector<double> maps_dx, vector<double> maps_dy,
                     double max_s)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1) % maps_x.size();

	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);
    double s_to_next_wp = maps_s[wp2] - s;
    if (s_to_next_wp < 0) {
        s_to_next_wp += max_s;
    }
    vector<double> normal_vector = interopolate_normal_vector({maps_dx[prev_wp], maps_dy[prev_wp]},
                                                              {maps_dx[wp2], maps_dy[wp2]},
                                                                s_to_next_wp, seg_s);
    double x = maps_x[prev_wp];
    double y = maps_y[prev_wp];
    double n_x = maps_x[wp2] - maps_x[prev_wp];
    double n_y = maps_y[wp2] - maps_y[prev_wp];
    double wp_distance = sqrt(n_x * n_x + n_y * n_y);
    x += seg_s * n_x / wp_distance;
    y += seg_s * n_y / wp_distance;

    x += d * normal_vector[0];
    y += d * normal_vector[1];

	return {x,y};

}

vector<double> predict_xy_position(double x, double y, double vx, double vy, double delta_t) {
    x += vx * delta_t;
    y += vy * delta_t;
    return {x, y};
}

vector<vector<double>> generate_trajectory(double car_x, double car_y, double car_yaw, double car_s,
                                           int lane,
                                           int old_lane,
                                           vector<double> previous_path_x,
                                           vector<double> previous_path_y,
                                           double step_time_interval, double ref_vel,
                                           vector<double> map_waypoints_x,
                                           vector<double> map_waypoints_y, vector<double> map_waypoints_s,
                                            vector<double> map_waypoints_dx, vector<double> map_waypoints_dy,
                                           double max_s) {
    vector<double> ptsx;
    vector<double> ptsy;

    // reference x, y, yaw states
    // either we will reference the starting point as where the car is or at the previous path and point
    double ref_x = car_x;
    double ref_y = car_y;
    double ref_yaw = deg2rad(car_yaw);

    // if previous size is almost empty, use the car as starting reference
    int prev_size = previous_path_x.size();
    if (prev_size < 2) {
        // Use two points that make the path tangent to the car
        double prev_car_x = car_x - cos(car_yaw);
        double prev_car_y = car_y - sin(car_yaw);

        ptsx.push_back(prev_car_x);
        ptsx.push_back(car_x);

        ptsy.push_back(prev_car_y);
        ptsy.push_back(car_y);
    } else { // use the previous path's end point as starting reference
        // Redefine reference state as previous path end point
        ref_x = previous_path_x[prev_size - 1];
        ref_y = previous_path_y[prev_size - 1];
        // ref_x = previous_path_x[1];
        // ref_y = previous_path_x[1];

        double ref_x_prev = previous_path_x[prev_size - 2];
        double ref_y_prev = previous_path_y[prev_size - 2];
        ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

        // Use two points that make the path tangent to the previous path's end point
        ptsx.push_back(ref_x_prev);
        ptsx.push_back(ref_x);

        ptsy.push_back(ref_y_prev);
        ptsy.push_back(ref_y);
    }

    // In frenet add evenly 30m spaced points ahead of the starting reference
    double planning_length = 30;
    /*if (lane != old_lane) {
        planning_length = 32;
    }*/
    vector<double> next_wp0 = getXY(car_s + planning_length, (2 + 4 * lane), map_waypoints_x, map_waypoints_y,
                                    map_waypoints_s,
                                    map_waypoints_dx, map_waypoints_dy, max_s);
    vector<double> next_wp1 = getXY(car_s + planning_length * 2, (2 + 4 * lane), map_waypoints_x, map_waypoints_y,
                                    map_waypoints_s, map_waypoints_dx, map_waypoints_dy, max_s);
    vector<double> next_wp2 = getXY(car_s + planning_length * 3, (2 + 4 * lane), map_waypoints_x, map_waypoints_y,
                                    map_waypoints_s, map_waypoints_dx, map_waypoints_dy, max_s);

    ptsx.push_back(next_wp0[0]);
    ptsx.push_back(next_wp1[0]);
    ptsx.push_back(next_wp2[0]);

    ptsy.push_back(next_wp0[1]);
    ptsy.push_back(next_wp1[1]);
    ptsy.push_back(next_wp2[1]);

    for (int i = 0; i < ptsx.size(); i++) {
        // shift car reference angle to 0 degrees
        double shift_x = ptsx[i] - ref_x;
        double shift_y = ptsy[i] - ref_y;

        ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
        ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
    }

    for (int i = 1; i < ptsx.size(); i++)
    {
        if (ptsx[i - 1] > ptsx[i])
        {
            double temp = ptsx[i];
            double temp_y = ptsy[i];

            int j = i;
            while (j > 0 && ptsx[j - 1] > temp)
            {
                ptsx[j] = ptsx[j - 1];
                ptsy[j] = ptsy[j - 1];
                j--;
            }
            ptsx[j] = temp;
            ptsy[j] = temp_y;
        }
    }

    // create a spline
    tk::spline s;
    for (int i = 0; i < ptsx.size(); i++) {
        if (ptsx[i] == ptsx[i + 1]) {
            ptsx[i + 1] += 0.1;
        }
    }
    // set (x,y) points to the spline
    // cout <<":" << ptsx[0] << " " << ptsx[1] << " " << ptsx[2] << endl;
    s.set_points(ptsx, ptsy);

    // Define the actual (x, y) points we will use for the planner




    vector<double> next_x_vals;
    vector<double> next_y_vals;

    for (int i = 0; i < previous_path_x.size(); i++) {
        next_x_vals.push_back(previous_path_x[i]);
        next_y_vals.push_back(previous_path_y[i]);
    }

    // Calculate how to break up spline points so that we travel at our desired reference velocity
    double target_x = 30.0;
    double target_y = s(target_x);
    double target_dist = sqrt(target_x * target_x + target_y * target_y);

    double x_add_on = 0;

    // Fill up the rest of our path planner after filling it with previous points, here we will always
    // output 50 points
    for (int i = 1; i <= 30 - previous_path_x.size(); i++) {
        double N = target_dist / (step_time_interval * ref_vel / 2.24);
        double x_point = x_add_on + target_x / N;
        double y_point = s(x_point);

        x_add_on = x_point;

        double x_ref = x_point;
        double y_ref = y_point;

        // rotate back to normal after rotating it earlier
        x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
        y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);

        x_point += ref_x;
        y_point += ref_y;

        next_x_vals.push_back(x_point);
        next_y_vals.push_back(y_point);

    }
    return {next_x_vals, next_y_vals};
}

vector<vector<vector<double>>> get_predictions(vector<vector<double>> sensor_fusion, int pred_steps,
                                               double step_time_interval) {
    vector<vector<vector<double>>> predictions;

    for (int i = 0; i < pred_steps; i++) {
        vector<vector<double>> predictions_one_step;
        for (int j = 0; j < sensor_fusion.size(); j++) {
            auto id = sensor_fusion[j][0];
            double x = sensor_fusion[j][1];
            double y = sensor_fusion[j][2];
            double vx = sensor_fusion[j][3];
            double vy = sensor_fusion[j][4];
            double v = sqrt(vx * vx + vy * vy);
            double s = sensor_fusion[j][5];
            double d = sensor_fusion[j][6];
            s += v * i * step_time_interval;

            // x += vx * i * step_time_interval;
            // y += vy * i * step_time_interval;
            predictions_one_step.push_back({id, s, d});

        }
        predictions.push_back(predictions_one_step);
    }
    return predictions;
}

vector<double> get_sd_prediction(vector<double> xy_prediction, vector<double> maps_x, vector<double> maps_y,
        vector<double> maps_s, vector<double> maps_dx, vector<double> maps_dy) {
    double x = xy_prediction.at(1);
    double y = xy_prediction.at(2);
    return getFrenet(x, y, maps_x, maps_y, maps_s, maps_dx, maps_dy);
}

bool check_collision(double s, double other_prev_s, double other_s) {
    // cout << "s: " << s << "\tother prev: " << other_prev_s << "\tother: " << other_s << endl;
    if (other_prev_s <= s - 20) {
        if (other_s >= s - 20) {
            return true;
        }
        return false;
    }
    if (other_prev_s > s + 20) {
        if (other_s <= s + 20) {
            return true;
        }
        return false;
    } else {
        return true;
    }

}
vector<double> get_helper_data(vector<vector<double>> trajectory, vector<vector<vector<double>>> predictions, int prev_size,
                     vector<double> maps_x, vector<double> maps_y, vector<double> maps_s, vector<double> maps_dx,
                     vector<double> maps_dy,
                     double step_time_interval) {
    // int PLANNING_HORIZON = 3;
    vector<double> ptsx = trajectory.at(0);
    vector<double> ptsy = trajectory.at(1);
    double current_x = ptsx.at(0);
    double current_y = ptsy.at(0);
    vector<double> current_sd = getFrenet(current_x, current_y, maps_x, maps_y, maps_s, maps_dx, maps_dy);
    double end_x = ptsx.at(ptsx.size() - 1);
    double end_y = ptsy.at(ptsx.size() - 1);
    vector<double> end_sd = getFrenet(end_x, end_y, maps_x, maps_y, maps_s, maps_dx, maps_dy);
    double avg_speed = (end_sd.at(0) - current_sd.at(0)) / (ptsx.size() * step_time_interval);
    double d_change = abs(end_sd.at(1) - current_sd.at(1));

    double closest_approach = 999999;
    bool collides = false;
    double collides_at = 9999999;
    double last_snap_x = ptsx[0];
    double last_snap_y = ptsy[0];
    for (int i = 2; i < ptsx.size(); i++ ) {
        double x = ptsx[i];
        double y = ptsy[i];
        vector<double> sd = getFrenet(x, y, maps_x , maps_y, maps_s, maps_dx, maps_dy);
        double s = sd.at(0);
        double d = sd.at(1);
        int num_vehicles = predictions.at(0).size();
        int index_diff = prev_size < 2 ? 0: prev_size;
        for (int j = 0; j < num_vehicles; j++) {

            vector<double> sd_state = predictions.at(i - 1 + index_diff).at(j);
            vector<double> prev_sd_state = predictions.at(i - 2 + index_diff).at(j);
           //  vector<double> sd_state = get_sd_prediction(state, maps_x, maps_y, maps_s, maps_dx, maps_dy);
            // vector<double> prev_sd_state = get_sd_prediction(prev_state, maps_x, maps_y, maps_s, maps_dx, maps_dy);
            if (sd_state.at(2) >= d - 4 && sd_state.at(2) <= d + 4) {
                bool vehicle_collides = check_collision(s, prev_sd_state.at(1), sd_state.at(1));
                if (vehicle_collides) {
                    collides = true;
                    collides_at = 0.0;
                }
                double dist = abs(sd_state.at(1) - s);
                if (dist < closest_approach) {
                    closest_approach = dist;
                }
            }

        }
    }
    return {collides_at, closest_approach, avg_speed, d_change};
}

double get_collision_cost(vector<double> data) {
    double collides_at = data.at(0);
    // collides_at = max(collides_at - 8.0, 0.0);
    double cost = exp(-collides_at);
    cout << "collides_at : " << collides_at << ", cost is: " << cost << endl;
    return cost;

}

double get_buffer_cost(vector<double> data, double DESIRED_BUFFER) {
    double closest_approach = data.at(1);
    if (closest_approach < 20) {
        return 10.0;
    }
    double timesteps_away = closest_approach / data.at(2);
    if (timesteps_away > DESIRED_BUFFER) {
        return 0.0;
    }

    double multiplier = 1.0 - pow((timesteps_away / DESIRED_BUFFER), 2);
    return multiplier;

}

double get_inefficiency_cost(vector<double> data, double target_speed) {

    double speed = data.at(2);

    double diff = target_speed - speed;
    double pct = diff / target_speed;
    double multiplier = pow(pct, 2);
    return multiplier;
}

double get_change_lane_cost(vector<double> data) {
    return pow(data.at(3), 3);
}

double calculate_cost(vector<vector<double>> trajectory, vector<vector<vector<double>>> predictions, int prev_size,
vector<double> maps_x, vector<double> maps_y, vector<double> maps_s, vector<double> maps_dx, vector<double> maps_dy,
                      double step_time_interval) {
    double COLLISION  = pow(10, 6);
    double DANGER     = pow(10, 5);
    double REACH_GOAL = pow(10, 5);
    double COMFORT    = pow(10, 4);
    double EFFICIENCY = pow(10, 2);

    double DESIRED_BUFFER = 3;
    double TARGET_SPEED = 49.5;
    vector<double> data = get_helper_data(trajectory, predictions, prev_size, maps_x, maps_y, maps_s, maps_dx,
    maps_dy, step_time_interval);
    double collision_cost = COLLISION * get_collision_cost(data);
    double danger_cost = DANGER * get_buffer_cost(data, DESIRED_BUFFER);
    double inefficiency_cost = EFFICIENCY * get_inefficiency_cost(data, TARGET_SPEED);
    double change_lane_cost = COMFORT * get_change_lane_cost(data);
    return collision_cost + danger_cost + inefficiency_cost + change_lane_cost;

}

vector<string> get_possible_state(string state, int lane, int keep_lane_counter) {
    if (keep_lane_counter < 10) {
        return {"KL"};
    }
    if (state == "KL") {
        if (lane == 0) {
            return {"KL", "LCR"};
        }
        if (lane == 2) {
            return {"KL", "LCL"};
        }
        else {

            return {"KL", "LCL", "LCR"};
        }
    }
    else {
        return {"KL"};
    }
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  double step_time_interval = 0.02;
    string state = "KL";
    vector<string> states = {"KL", "LCL", "LCR"};
    int keep_lane_counter = 0;


  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
    int lane = 1;
    double ref_vel = 0.0;

  h.onMessage([&lane, &ref_vel, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,
  &step_time_interval, &max_s, &state, &keep_lane_counter](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data

          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            int prev_size = previous_path_x.size();
            // cout << "prev size is: " << prev_size << endl;
            if (prev_size > 0) {
                car_s = end_path_s;
            }

            // prediction
            // prediction steps
            int pred_steps = 60;
            vector<vector<vector<double>>> predictions = get_predictions(sensor_fusion, pred_steps, step_time_interval);


            int left_lane = lane - 1;
            int right_lane = lane + 1;
            bool too_close = false;
            bool can_change_left = false;
            bool can_change_right = false;
            if (left_lane >= 0) {
                can_change_left = true;
            }
            if (right_lane <= 2) {
                can_change_right = true;
            }

            // find ref_v to use
            // vector<vector<double>> vector_snap = predictions[prev_size - 1];

            for (int i = 0; i < sensor_fusion.size(); i++) {

                double s = predictions.at(prev_size).at(i).at(1);
                double d = predictions.at(prev_size).at(i).at(2);
                // vector<double> sd_position = getFrenet(x, y, map_waypoints_x, map_waypoints_y,
                // map_waypoints_s, map_waypoints_dx, map_waypoints_dy);
                // double s = sd_position[0];
                // float d = sd_position[1];
                if (d < (2 + 4 * lane + 2) && d > (2 + 4 * lane - 2)) {
                    if ((s > car_s) && (s - car_s) < 30) {
                        too_close = true;

                    }
                }
                if (can_change_left) {

                    if ( d < (2 + 4 * left_lane + 2) && d > (2 + 4 * left_lane - 2)) {
                        if ((s > car_s) && (s - car_s) < 30) {
                            can_change_left = false;
                        }
                        if ((s < car_s) && (car_s - s) < 30) {
                            can_change_left = false;
                        }
                    }
                }
                if (can_change_right) {
                    if (d < (2 + 4 * right_lane + 2) && d > (2 + 4 * right_lane - 2)) {
                        if ((s > car_s) && (s - car_s) < 30) {
                            can_change_right = false;
                        }
                        if ((s < car_s) && (car_s - s) < 30) {
                            can_change_right   = false;
                        }
                    }
                }

            }
            if (too_close) {
                if (can_change_left) {
                    lane -= 1;
                } else if (can_change_right){
                    lane += 1;
                } else {
                    ref_vel -= 0.224;
                }
            } else {
                if (ref_vel < 49.5) {
                    ref_vel += 0.224;
                }
            }

            auto best_trajectory = generate_trajectory(car_x, car_y, car_yaw, car_s,
                                                       lane,
                                                       lane,
                                                       previous_path_x,
                                                       previous_path_y,
                                                       step_time_interval,
                                                       ref_vel,
                                                       map_waypoints_x, map_waypoints_y,
                                                       map_waypoints_s,
                                                       map_waypoints_dx, map_waypoints_dy, max_s);


            //  double min_cost = 99999999999999999; // very big number
            // vector<vector<double>> best_trajectory;
            /* if (too_close) {
                // slow down or change lane
                vector<string> possible_state = get_possible_state(state, lane, keep_lane_counter);
                //vector<int> possible_lanes = get_possible_lanes(lane);
                int old_lane = lane;

                for (int i = 0; i < possible_state.size(); i++) {
                    string new_state = possible_state[i];
                    int proposed_lane = lane;
                    if (new_state == "LCL") {
                        proposed_lane = lane - 1;
                    } else if (new_state == "LCR") {
                        proposed_lane = lane + 1;
                    }
                    /*if (proposed_lane == lane) {
                        ref_vel -= 0.224;
                    }
                    vector<vector<double>> trajectory = generate_trajectory(car_x, car_y, car_yaw, car_s,
                                                                            proposed_lane,
                                                                            lane,
                                                                            previous_path_x,
                                                                            previous_path_y,
                                                                            step_time_interval,
                                                                            ref_vel,
                                                                            map_waypoints_x, map_waypoints_y,
                                                                            map_waypoints_s,
                                                                            map_waypoints_dx, map_waypoints_dy, max_s);

                    double cost = calculate_cost(trajectory, predictions, prev_size, map_waypoints_x,
                    map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy, step_time_interval);
                    if (cost < min_cost) {
                        min_cost = cost;
                        best_trajectory = trajectory;
                        lane = proposed_lane;
                        state = new_state;

                    }
                }
                if (min_cost >= pow(10, 6)) { // all solutions danger, just slow down
                    ref_vel -= 0.224;

                    state = "KL";
                    lane = old_lane;
                    best_trajectory = generate_trajectory(car_x, car_y, car_yaw, car_s,
                                                          old_lane,
                                                          old_lane,
                                                          previous_path_x,
                                                          previous_path_y,
                                                          step_time_interval,
                                                          ref_vel,
                                                          map_waypoints_x, map_waypoints_y,
                                                          map_waypoints_s,
                                                          map_waypoints_dx, map_waypoints_dy, max_s);

                }
                if (state == "KL") {
                    keep_lane_counter++;
                } else {
                    keep_lane_counter = 0;
                }


                // ref_vel -= 0.224;
            } else {
                if (ref_vel < 49.5) {
                    ref_vel += 0.224;
                }
                state = "KL";
                keep_lane_counter++;
                vector<vector<double>> trajectory = generate_trajectory(car_x, car_y, car_yaw, car_s, lane, lane,
                                                                    previous_path_x,
                                                                    previous_path_y,
                                                                    step_time_interval,
                                                                    ref_vel,
                                                                    map_waypoints_x, map_waypoints_y,
                                                                    map_waypoints_s,
                                                                    map_waypoints_dx, map_waypoints_dy, max_s);
                best_trajectory = trajectory;

            } */



            // Create a list of widely spaced (x, y) waypoints, evenly spaced at 30m
            // Later we will interoplate these waypoints with a spline and fill it in with more points


            vector<double> next_x_vals = best_trajectory.at(0);
            vector<double> next_y_vals = best_trajectory.at(1);
            json msgJson;
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

















































































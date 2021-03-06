/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>


#include "particle_filter.h"

#define NUMBER_OF_PARTICLES 500 // Can be decreased (even 12 particles can pass the test)
#define EPS 0.001  // Just a small number

using namespace std;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    particles.resize(num_particles);
    weights.resize(num_particles);
    
    default_random_engine gen;
    normal_distribution<double> dist_x (x,std[0]);
    normal_distribution<double> dist_y (y,std[1]);
    normal_distribution<double> dist_theta (theta,std[2]);
    
    for (int i=0; i < particles.size(); i++){
        
        particles[i].id = i;
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
        particles[i].weight = 1;
        weights[i] = particles[i].weight;
        
    }
    
    is_initialized = true;
    
    
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    
    
    default_random_engine gen;
    normal_distribution<double> dist_x (0,std_pos[0]);
    normal_distribution<double> dist_y (0,std_pos[1]);
    normal_distribution<double> dist_theta (0,std_pos[2]);
    
    const double yawr_dt = yaw_rate*delta_t;
    
    
    for (int i=0; i<particles.size(); i++){
        
        //Prediction Step
        if (fabs(yaw_rate) > 0.001){
            const double vel_yawr = velocity/yaw_rate;
            particles[i].x += vel_yawr*(sin(particles[i].theta + yawr_dt)-sin(particles[i].theta));
            particles[i].y += vel_yawr*(cos(particles[i].theta) - cos(particles[i].theta + yawr_dt));
            particles[i].theta +=  yawr_dt;
            
        }else{
            
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);
            
            
            
        }
        
        //Adding Noise
        particles[i].x = particles[i].x + dist_x(gen);
        particles[i].y = particles[i].y + dist_y(gen);
        particles[i].theta = particles[i].theta + dist_theta(gen);
        
        
    }
    
}

inline void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.
    
    for (int i=0; i < observations.size(); i++) {
        
        std::vector<double> distance;
        distance.resize(predicted.size());
        std::vector<int> obs_id;
        obs_id.resize(predicted.size());
        double obs_x = observations[i].x;
        double obs_y = observations[i].y;
        
        
        for (int j=0; j < predicted.size(); j++) {
            
            double pre_x = predicted[j].x;
            double pre_y = predicted[j].y;
            
            distance[j] = dist(obs_x, obs_y, pre_x, pre_y);
            obs_id[j] = predicted[j].id;
            
        }//end j
        
        
        auto min_value = min_element(begin(distance), end(distance));\
        for (int id=0; id < distance.size() ; id++){ //Assign Id to observations
            if (*min_value == distance[id]){
                observations[i].id = obs_id[id];
                id = distance.size();
            }
        }
        
        
        
        
    }//end i
    
    
    
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html
    
    std::vector<LandmarkObs> car_to_map_obs;
    car_to_map_obs.resize(observations.size());
    std::vector <double>  distance ; // distance between observation and landmark
    distance.resize(map_landmarks.landmark_list.size());
    
    const double gauss_norm = 1.0/(2*M_PI*std_landmark[0]*std_landmark[1]);
    const double sigma_xx = 2*pow(std_landmark[0],2);
    const double sigma_yy = 2*pow(std_landmark[1],2);
    double sum_weights = 0.0;
    
    
    
    
    
    for (int z = 0 ; z < num_particles ; z++){
        
        
        double par_x = particles[z].x;
        double par_y = particles[z].y;
        
        std::vector<LandmarkObs> predicted;
        
        for (int i=0 ; i < map_landmarks.landmark_list.size(); i++){ //Particles Prediction in range
            double lm_x = map_landmarks.landmark_list[i].x_f;
            double lm_y = map_landmarks.landmark_list[i].y_f;
            int lm_id = map_landmarks.landmark_list[i].id_i;
            
            if (dist(par_x,par_y,lm_x,lm_y) <= sensor_range){
                
                predicted.push_back(LandmarkObs{lm_id,lm_x,lm_y});
                
                
            }//end if
            
        }//end i Particles Prediction in sensor range
        
        
        for (int j=0; j < observations.size() ;j++){ //Car to Map Coordinates
            const double cos_ = cos(particles[z].theta);
            const double sin_ = sin(particles[z].theta);
                                    
            
            car_to_map_obs[j].x = par_x + (cos_*observations[j].x) - (sin_*observations[j].y);
            car_to_map_obs[j].y = par_y + (sin_*observations[j].x) + (cos_*observations[j].y);
            
        }//end j Car to Map coordinates*/
        
        
        dataAssociation(predicted, car_to_map_obs);
        
        //Set Association Test
        std::vector<int> associations_id;
        std::vector<double> sense_x;
        std::vector<double> sense_y;
        
        
        for (int t=0; t < car_to_map_obs.size(); t++){
            associations_id.push_back(car_to_map_obs[t].id);
            sense_x.push_back(car_to_map_obs[t].x);
            sense_y.push_back(car_to_map_obs[t].y);
        }
        
        SetAssociations(particles[z], associations_id, sense_x, sense_y);
        
        double obs_x=0.0;
        double obs_y=0.0;
        double land_x=0.0;
        double land_y=0.0;
        
        double eq_term1 = 0.0 ;
        double eq_term2 = 0.0 ;
        long double obs_weight = 0.0 ;
        particles[z].weight  = 1.0;
        
        
        for (int j = 0; j < car_to_map_obs.size() ; j++){
            
            obs_x = car_to_map_obs[j].x;
            obs_y = car_to_map_obs[j].y;
            
            
            for(int i = 0; i < map_landmarks.landmark_list.size(); i++){
                
                if (car_to_map_obs[j].id == map_landmarks.landmark_list[i].id_i){
                    land_x = map_landmarks.landmark_list[i].x_f;
                    land_y = map_landmarks.landmark_list[i].y_f;
                }//end if
            }//end i landmarks
            
            /*obs_x = 0;
             obs_y = 5;
             land_x = 2;
             land_y  = 1;*/
            
            eq_term1 = (pow(obs_x-land_x,2))/(sigma_xx);
            eq_term2 = (pow(obs_y-land_y,2))/(sigma_yy);
            obs_weight = gauss_norm*exp(-1.0*(eq_term1+eq_term2));
    
            
            particles[z].weight *= obs_weight;
                
            
            
            /*
            cout << "======================================" << endl;
            cout << "OBSERVATION " << j << endl;
            cout << "Obs X: " << obs_x << " Obs Y: " << obs_y << endl;
            cout << "Land_X: " << land_x << " Land_Y: " << land_y << endl;
            cout << "Weights: " << obs_weight << endl;
            */
            
            
        }// end j car_to_maps_obs
        
        obs_weight = 0.0;
        weights[z] = particles[z].weight;
        sum_weights += weights[z];
        
        
        
    } //end for z main loop
    
    /*
    for (int i = 0; i < num_particles; i++){ //Regularization
        particles[i].weight /= sum_weights;
        weights[i] = particles[i].weight;
        
        
    }*/
    
    sum_weights = 0.0;
    
    
    
} // End of Update Function

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    vector<Particle> new_particles (num_particles);
    random_device rd;
    default_random_engine gen(rd());
    discrete_distribution<int> index(weights.begin(), weights.end());
    
    for (int i = 0; i < num_particles; ++i) {
        
        new_particles[i] = particles[index(gen)];
        
    }
    
    particles = new_particles;
    
    
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                         const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
    
    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();
    
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    
    return particle;
    
    
    
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
    vector<double> v = best.sense_x;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
    vector<double> v = best.sense_y;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

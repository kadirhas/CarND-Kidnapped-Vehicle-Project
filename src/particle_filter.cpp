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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 100;
	double wtemp = 1.0;
	default_random_engine gen;
	normal_distribution<double> dist_x(x,std[0]);
	normal_distribution<double> dist_y(y,std[1]);
	normal_distribution<double> dist_theta(theta,std[2]);
	for (int i = 0; i < num_particles; i++)
	{
		Particle temp;
		//std[] = sigma_pos [3] = {0.3, 0.3, 0.01}; // GPS measurement uncertainty [x [m], y [m], theta [rad]]
		temp.x = dist_x(gen);
		temp.y = dist_y(gen);
		temp.theta = dist_theta(gen);
		temp.weight = wtemp;
		temp.id = i;
		particles.push_back(temp);
		weights.push_back(wtemp);
	}
	is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	//std[] = sigma_pos [3] = {0.3, 0.3, 0.01}; // GPS measurement uncertainty [x [m], y [m], theta [rad]]
	default_random_engine gen;
	normal_distribution<double> dist_x(0,std_pos[0]);
	normal_distribution<double> dist_y(0,std_pos[1]);
	normal_distribution<double> dist_theta(0,std_pos[2]);
	if (fabs(yaw_rate)<0.0001)
	{
		for (int i = 0; i<num_particles; i++)
		{
			//you also need to add std deviation
			particles[i].x += velocity*delta_t*cos(particles[i].theta);
			particles[i].y += velocity*delta_t*sin(particles[i].theta);
		}
	}
	else
	{
		const double yawdt = yaw_rate*delta_t;
		for (int i = 0; i<num_particles; i++)
		{
			double theta_0 = particles[i].theta;
			//add noise
			particles[i].x += (velocity/yaw_rate)*(sin(theta_0+yawdt) - sin(theta_0));
			particles[i].y += (velocity/yaw_rate)*(cos(theta_0) - cos(theta_0+yawdt));
			particles[i].theta += yawdt;
		}
	}
	for (int i = 0; i<num_particles; i++)
	{
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);
	}

}
// a
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int obv = 0; obv < observations.size(); obv++)
	{
		double min_dist = 10000000000000000;
		int prdid = -1;
		for (int pr = 0; pr<predicted.size(); pr ++)
		{
			double distance = dist(observations[obv].x, observations[obv].y, predicted[pr].x, predicted[pr].y);
			if (distance<min_dist)
			{
				min_dist = distance;
				prdid = predicted[pr].id;
			}
		}

		observations[obv].id = prdid;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
		//std_landmark = sigma_landmark [2] = {0.3, 0.3}; // Landmark measurement uncertainty [x [m], y [m]]
		for (int i = 0; i < num_particles; i++)
		{
			double px = particles[i].x;
			double py = particles[i].y;
			double ptheta = particles[i].theta;
			//initiate particle weight for calculation
			particles[i].weight = 1.0;
			vector<LandmarkObs> obsPredictions;
			// only consider the landmarks in the sensor range
			for (int j = 0; j<map_landmarks.landmark_list.size(); j++)
			{
				float lmx = map_landmarks.landmark_list[j].x_f;
				float lmy = map_landmarks.landmark_list[j].y_f;
				if (dist(px,py,lmx,lmy)<sensor_range)
				{
					obsPredictions.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i,lmx,lmy});
				}
			}
			//transform sensor obstacles to global coordinates
			vector<LandmarkObs> obsGlobal;

			for (int obs = 0; obs < observations.size(); obs++ )
			{
				LandmarkObs temp;
				temp.x = observations[obs].x * cos(ptheta) - observations[obs].y * sin(ptheta) + px;
				temp.y = observations[obs].x * sin(ptheta) + observations[obs].y * cos(ptheta) + py;
				temp.id = observations[obs].id;
				obsGlobal.push_back(temp);
			}
			//associate them for each particle
			dataAssociation(obsPredictions,obsGlobal);
			//calculate the weights for each landmark
			for (int obs = 0; obs < obsGlobal.size(); obs++)
			{
				double predictionX;
				double predictionY;
				int currentID = obsGlobal[obs].id;
				//get the matched prediction
				for (int match = 0; match < obsPredictions.size(); match++)
				{
					if (currentID==obsPredictions[match].id)
					{
						predictionX = obsPredictions[match].x;
						predictionY = obsPredictions[match].y;
						break;
					}
				}//for (int match = 0; match < obsPredictions.size(); match++)
				//weight calculation and production with the earlier ones
				particles[i].weight *= exp( -(pow(predictionX-obsGlobal[obs].x,2)/(2*pow(std_landmark[0],2))+(pow(predictionY-obsGlobal[obs].y,2)/(2*pow(std_landmark[1],2)) )) )/(2*M_PI*std_landmark[0]*std_landmark[1]);
			}//for (int obs = 0; obs < obsGlobal.size(); obs++)
			weights[i] = particles[i].weight;
		} //end of particle loop
	} // end of update weights function

	void ParticleFilter::resample() {
		// TODO: Resample particles with replacement with probability proportional to their weight.
		// NOTE: You may find std::discrete_distribution helpful here.
		//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
		vector<Particle> newParticles;
		uniform_int_distribution<int> newIndex(0,num_particles-1);
		uniform_real_distribution<double> dist(0, 1);
		default_random_engine gen;
		int index = newIndex(gen);
		double beta = 0.0;
		double maxWeight = *max_element(weights.begin(), weights.end());
		for (int i = 0; i<num_particles;i++)
		{
			beta += dist(gen) * 2.0 * maxWeight;
			while(beta>weights[index])
			{
				beta -= weights[index];
				index = (index + 1)%num_particles;
			}
			newParticles.push_back(particles[index]);
		}
		particles = newParticles;
	}

	Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
		const std::vector<double>& sense_x, const std::vector<double>& sense_y)
		{
			//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
			// associations: The landmark id that goes along with each listed association
			// sense_x: the associations x mapping already converted to world coordinates
			// sense_y: the associations y mapping already converted to world coordinates

			particle.associations= associations;
			particle.sense_x = sense_x;
			particle.sense_y = sense_y;
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

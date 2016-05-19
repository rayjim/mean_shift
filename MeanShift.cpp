#include <stdio.h>
#include <math.h>
#include "MeanShift.h"

#define EPSILON 0.0000001
#define GROUP_DISTANCE_TOLERANCE 0.1 
#define MAXIMUM_DISTANCE 100000000

double euclidean_distance(const vector<double> &point_a, const vector<double> &point_b){
    double total = 0;
    for(int i=0; i<point_a.size(); i++){
        total += (point_a[i] - point_b[i]) * (point_a[i] - point_b[i]);
    }
    return sqrt(total);
}

double gaussian_kernel(double distance, double kernel_bandwidth){
    double temp =  exp(-(distance*distance) / (kernel_bandwidth));
    return temp;
}

// To grouping the final results of the meanshift
// input: the results of the shifted points
// output: A vector of the labels of the shifted points

vector<int > group_points(vector<vector<double> > points){
    vector<vector<double> > temp_points; //save the points for the peak
    vector<int> indx_points;
    //search for all points
    for (int i = 0; i < points.size();i++){
        //compared points with points in the temp list
	double min_distance = GROUP_DISTANCE_TOLERANCE;
        int min_index = -1;	
	for (int j = 0; j < temp_points.size(); j++){
	    double dist = euclidean_distance(points[i], temp_points[j]);
	    if (dist < min_distance){
		    min_distance = dist;
	    	    min_index = j;
	    }
	}
	//if depends on the index,
	if (min_index != -1){
          indx_points.push_back(min_index);  
	}
	else{//added to the new temp_points
	   indx_points.push_back(temp_points.size()+1);
           temp_points.push_back(points[i]);
	}
    }	
    return indx_points;
}

void MeanShift::set_kernel( double (*_kernel_func)(double,double) ) {
    if(!_kernel_func){
        kernel_func = gaussian_kernel;
    } else {
        kernel_func = _kernel_func;    
    }
}

vector<double> MeanShift::shift_point(const vector<double> &point, const vector<vector<double> > &points, double kernel_bandwidth) {
    vector<double> shifted_point = point;
    for(int dim = 0; dim<shifted_point.size(); dim++){
        shifted_point[dim] = 0;
    }
    double total_weight = 0;
    for(int i=0; i<points.size(); i++){
        vector<double> temp_point = points[i];
        double distance = euclidean_distance(point, temp_point);
        double weight = kernel_func(distance, kernel_bandwidth);
        for(int j=0; j<shifted_point.size(); j++){
            shifted_point[j] += temp_point[j] * weight;
        }
        total_weight += weight;
    }

    for(int i=0; i<shifted_point.size(); i++){
        shifted_point[i] /= total_weight;
    }
    return shifted_point;
}

vector<int > MeanShift::cluster(vector<vector<double> > points, double kernel_bandwidth){
    vector<bool> stop_moving;
    stop_moving.reserve(points.size());
    vector<vector<double> > shifted_points = points;
    double max_shift_distance;
    do {
        max_shift_distance = 0;
        for(int i=0; i<shifted_points.size(); i++){
            if (!stop_moving[i]) {
                vector<double>point_new = shift_point(shifted_points[i], points, kernel_bandwidth);
                double shift_distance = euclidean_distance(point_new, shifted_points[i]);
                if(shift_distance > max_shift_distance){
                    max_shift_distance = shift_distance;
                }
                if(shift_distance <= EPSILON) {
                    stop_moving[i] = true;
                }
                shifted_points[i] = point_new;
            }
        }
        printf("max_shift_distance: %f\n", max_shift_distance);
    } while (max_shift_distance > EPSILON);
    vector<int> labels = group_points(shifted_points);
    return labels;
}


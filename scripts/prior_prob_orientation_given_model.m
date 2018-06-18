function prob = prior_prob_orientation_given_model(...
    estimated_orientation, given_orientation, max_angle_err, max_shift_err)
    
    esti_theta = estimated_orientation.theta;
    given_theta = given_orientation.theta;
    esti_shift = estimated_orientation.shift;
    given_shift = given_orientation.shift;
    
    prob = 0;
    
    if ((esti_theta - max_angle_err) <= given_theta) && (given_theta <= (esti_theta + max_angle_err))
        if ((esti_shift - max_shift_err) <= given_shift) && (given_shift <= (esti_shift + max_shift_err))
            prob = (1/(2*max_angle_err + 1))*(1/(2*max_shift_err + 1));
        end
    end
end
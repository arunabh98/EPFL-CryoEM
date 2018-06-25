function prob = prior_prob_orientation_given_model(...
    first_orientation, given_orientation, max_angle_err, max_shift_err)
    
    first_theta = first_orientation.theta;
    given_theta = given_orientation.theta;
    first_shift = first_orientation.shift;
    given_shift = given_orientation.shift;

    prob = 0;
    
    if ((first_theta - max_angle_err) <= given_theta) && (given_theta <= (first_theta + max_angle_err))
        if ((first_shift - max_shift_err) <= given_shift) && (given_shift <= (first_shift + max_shift_err))
            prob = (1/(2*max_angle_err + 1))*(1/(2*max_shift_err + 1));
        end
    end
end
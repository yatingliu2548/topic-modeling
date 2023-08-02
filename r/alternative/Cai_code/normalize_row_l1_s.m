function [ theta_l1, d_theta ] = normalize_row_l1_s( theta )
    d_theta = sum(theta,2);
    S = find(d_theta>0);
    theta_l1 = theta;
    theta_l1(S,:) = bsxfun(@times, theta(S,:), 1./(sum(theta(S,:), 2)));
end

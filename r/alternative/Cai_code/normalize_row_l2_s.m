function [ theta_l2 ] = normalize_row_l2_( theta )
    d_theta = sqrt(sum(theta'.*theta'));
    S = find(d_theta>0);
    theta_l2 = theta;
    theta_l2(S,:) = bsxfun(@times, theta(S,:), 1./(sqrt(sum(theta(S,:)'.*theta(S,:)')))');
end
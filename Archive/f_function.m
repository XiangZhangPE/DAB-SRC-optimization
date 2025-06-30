% Define f(Dy1, Dy2, Dp) function
function f = f_function(vars, r)
    Dy1 = vars(1);
    Dy2 = vars(2);
    Dp = vars(3);

    region_flag((Dy1>Dy2) & (Dp<(Dy1-Dy2)/2)) = 1;  % 区域 A
    region_flag((Dy2>=Dy1) & (Dp<(Dy2-Dy1)/2)) = 2;  % 区域 B
    region_flag((Dp>=abs(Dy1-Dy2)/2) & (Dp<(Dy2+Dy1)/2) & (Dp<(1-(Dy2+Dy1)/2))) = 3;  % 区域 C
    region_flag((((Dy2+Dy1)/2)<= Dp) & (Dp<(1-(Dy2+Dy1)/2))) = 4;  % 区域 D
    region_flag(((1-(Dy2+Dy1)/2)<= Dp) & (((Dy2+Dy1)/2)>= Dp)) = 5;  % 区域 E
    % region_flag(((1-(Dy2-Dy1)/2)<=Dp) & (Dy2>Dy1)) = 6;  % 区域 F
    % region_flag(((1-(Dy1-Dy2)/2)<=Dp) & (Dy2<=Dy1)) = 7;  % 区域 G
    % region_flag((Dp>=(1-(Dy1+Dy2)/2)) & (Dp>=(Dy1+Dy2)/2) & (Dp<(1-abs(Dy1-Dy2)/2))) = 8;  % 区域 H

    % Calculate f value
    switch region_flag
        case 1
            % Region A
            f = ((sin(pi*r.*Dy2/2).*sin(pi*r.*Dp).*cos(pi*r.*(1-Dy1)/2))./(2*r.*cos(pi*r/2)));
        case 2
            % Region B
            f = ((sin(pi*r.*Dy1/2).*sin(pi*r.*Dp).*cos(pi*r.*(1-Dy2)/2))./(2*r.*cos(pi*r/2)));
        case 3
            % Region C
            f = 0;
        case 4
            % Region D
            f = 0;
        case 5
            % Region E
            f = ((cos(pi*r*(Dp-0.5))*cos(pi*r*(1-Dy1)/2)*cos(pi*r*(1-Dy2)/2))/(pi*r*cos(pi*r/2))-1/(pi*r));
        % case 6
        %     % Region F
        %     g = ... % Add the specific formula for region F
        % case 7
        %     % Region G
        %     g = ... % Add the specific formula for region G
        % case 8
        %     % Region H
        %     g = ... % Add the specific formula for region H
        otherwise
            error('Unknown region');
    end

end
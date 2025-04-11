clc; clear; 
%%
xMesh = linspace(0, 10, 50);
Pr = 100;
global const;
const = 500;

% Solve f(x)
sol_f = solve_f(xMesh);
x_vals = sol_f.x;
f_vals = sol_f.y(1,:);
df_vals = sol_f.y(2,:);

% Solve theta 
solinit = bvpinit(x_vals, [0, 0]);
sol_theta = bvp4c(@(x,y) ode_dthetadx(x, y, sol_f, Pr), @bc_dthetadx, solinit);

figure(2);
plot(sol_theta.x, sol_theta.y(1,:), 'LineWidth', 2)
xlabel('x')
ylabel('\theta(x)')
title('Solution of Thermal Boundary Layer')
grid on
hold on


function sol_f = solve_f(xMesh)
    solinit = bvpinit(xMesh, [0, 0, 0.0]);  % y1=f, y2=f', y3=f''
    sol_f = bvp4c(@ode_dfdx, @bc_dfdx, solinit);

    x = sol_f.x;
    y = sol_f.y;
    figure(1)
    plot(x, y(1,:), '-', x, y(2,:), '--', x, y(3,:), '-.')
    legend('f(x)', 'df/dx', 'd^2f/dx^2')
    xlabel('x')
    ylabel('Solution')
    title('Blasius-like Equation Solution')
    grid on
    hold on
end

function dfdx = ode_dfdx(x, y)
    dfdx = zeros(3,1);
    dfdx(1) = y(2);                     % y1' = y2
    dfdx(2) = y(3);                     % y2' = y3
    dfdx(3) = -0.5 * y(1) * y(3);       % y3' = -0.5*f*f''
end

function out = bc_dfdx(ya, yb)
    out = [ya(1);         % f(0) = 0
           ya(2);         % f'(0) = 0
           yb(2) - 1];    % f'(∞) = 1
end

function dthetadx = ode_dthetadx(x, y, sol_f, Pr)
    f = deval(sol_f, x);     % [f; f'; f'']
    fx = f(1);
    dfdx = f(2);

    dthetadx = zeros(2,1);
    dthetadx(1) = y(2);
    dthetadx(2) = (-Pr/2) * (fx * y(2) - dfdx * y(1)); 
end

function out = bc_dthetadx(ya, yb)
    global const;
    out = [ya(2) + const;         % f'(0) = const
           yb(2)];            % f'(∞) = 0
end
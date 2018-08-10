function sx_checkFlight(filename)

% usage ex:
% sx_checkFlight('sea003.299.gli.sub.80')

[meta, data] = sx2mat(filename);

% plot angle
figure
plot(data(:,1)-data(1,1), data(:,5))
hold on
plot(data(:,1)-data(1,1), data(:,6), 'r')
ylabel('pitch/roll')
xlabel('Time (seconds in profiles)')
legend('pitch', 'roll')
set(gca, 'ygrid', 'on')
title([filename])


% plot velocity
figure
plot(data(2:end,1)-data(1,1), diff(data(:,7))./diff(data(:,1)))
ylabel('vert. vel.')
xlabel('Time (seconds in profiles)')
set(gca, 'ygrid', 'on')
title([filename])


% plot Depth
figure
plot(data(:,1)-data(1,1), data(:,7))
ylabel('Depth')
xlabel('Time (seconds in profiles)')
set(gca, 'ygrid', 'on')
set(gca, 'xgrid', 'on')
title([filename])

% plot Ballast
figure
plot(data(:,1)-data(1,1), data(:,13))
hold on
plot(data(:,1)-data(1,1), data(:,14), 'r')
legend('Cmd', 'Pos')
ylabel('Ballast')
xlabel('Time (seconds in profiles)')
set(gca, 'ygrid', 'on')
title([filename])

% plot Heading
figure
plot(data(:,1)-data(1,1), data(:,12))
hold on
plot(data(:,1)-data(1,1), data(:,4), 'r')
legend('Cmd', 'Pos')
ylabel('Heading')
xlabel('Time (seconds in profiles)')
set(gca, 'ygrid', 'on')
title([filename])

R = input('want more plots (NavState, Alarm, etc.)? (y/n)','s');
if strcmp(R, 'y')~=1
    return
end
    
% plot NavState
figure
plot(data(:,1)-data(1,1), data(:,2))
ylabel('NavState')
xlabel('Time (seconds in profiles)')
set(gca, 'ygrid', 'on')
title([filename])

% plot alarms
figure
plot(data(:,1)-data(1,1), data(:,3))
ylabel('Alarm')
xlabel('Time (seconds in profiles)')
set(gca, 'ygrid', 'on')
title([filename])
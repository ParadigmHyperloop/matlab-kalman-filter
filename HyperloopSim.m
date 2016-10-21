% Hyperloop pod simulation for the OpenLoop Team
% Created by Alex Goldstein and Adam Shaw

% function [] = HyperloopSim()
    disp('Simulation Started')
    
        
    globals = struct(globalData());
    pod = struct(podData());
    tube = struct(tubeData());
    
    failed=false;
    sensorType=[];
    sensorNum=[];
    sensorFailState=[];
    
    
    if globals.correctCovariance
         globals.peTopSIMCovConst = globals.peTopCovConst;
         globals.peLeftSIMCovConst = globals.peLeftCovConst;
         globals.peRightSIMCovConst = globals.peRightCovConst;
         globals.pitotSIMCovConst = globals.pitotCovConst;
         globals.distDownSIMCovConst = globals.distDownCovConst;
         globals.distDownRailSIMCovConst = globals.distDownRailCovConst;
         globals.distSideSIMCovConst = globals.distSideCovConst;
         globals.IMUAccelSIMCovConst = globals.IMUAccelCovConst;
         globals.IMUGyroSIMCovConst = globals.IMUGyroCovConst;
        
         globals.peTopSIMCovLin = globals.peTopCovLin;
         globals.peLeftSIMCovLin = globals.peLeftCovLin;
         globals.peRightSIMCovLin = globals.peRightCovLin;
         globals.pitotSIMCovLin = globals.pitotCovLin;
         globals.distDownSIMCovLin = globals.distDownCovLin;
         globals.distDownRailSIMCovLin = globals.distDownRailCovLin;
         globals.distSideSIMCovLin = globals.distSideCovLin;
         globals.IMUAccelSIMCovLin = globals.IMUAccelCovLin;
         globals.IMUGyroSIMCovLin = globals.IMUGyroCovLin;
               
         globals.peTopSIMCovZero = globals.peTopCovZero;
         globals.peLeftSIMCovZero = globals.peLeftCovZero;
         globals.peRightSIMCovZero = globals.peRightCovZero;
         globals.pitotSIMCovZero = globals.pitotCovZero;
         globals.distDownSIMCovZero = globals.distDownCovZero;
         globals.distDownRailSIMCovZero = globals.distDownRailCovZero;
         globals.distSideSIMCovZero = globals.distSideCovZero;
         globals.IMUAccelSIMCovZero = globals.IMUAccelCovZero;
         globals.IMUGyroSIMCovZero = globals.IMUGyroCovZero;
    end
    
    kalmanFreq = globals.kalmanTimestep/globals.timestep;

    % temporary vectors for debugging purposes
    rolls = zeros(1,globals.numSteps/kalmanFreq+1);
    sensorInfo=zeros(25,globals.numSteps/kalmanFreq+1);
    sensorRecord=zeros(25,globals.numSteps/kalmanFreq+1);
    scannerDistances = zeros(globals.numSteps/kalmanFreq+1);
    
    photoelectricCount=zeros(9,1);
    lastPE=-5*ones(9,1);
    prevPE=zeros(9,3);
    crossIn=zeros(9,1);

    
    peRecord = nan(9,globals.numSteps/kalmanFreq+1);
    noisyKalmanPosIMU = zeros(3,globals.numSteps/kalmanFreq+1);
    noisyKalmanVelIMU = zeros(3,globals.numSteps/kalmanFreq+1);
    noisyKalmanPosIMU(3,1) = .003;
    noisyKalmanRotIMU=[1 0 0; 0 1 0; 0 0 1];
    noisyKalmanQIMU=zeros(4,globals.numSteps/kalmanFreq+1);
    noisyKalmanQIMU(:,1)=[0.000001;0.000001;0.000001; sqrt(1-3*.000001^2)];
    gyroData= zeros(3,globals.numSteps/kalmanFreq+1);
    noisyGyroData= zeros(3,globals.numSteps/kalmanFreq+1);
    qdotIMU= zeros(4,globals.numSteps/kalmanFreq+1);
    % pod state variables
    transPos = zeros(3, globals.numSteps);
    transVel = zeros(3, globals.numSteps);
    transAcc = zeros(3, globals.numSteps);
    transPos(3,1) = .003;

    rotPos = zeros(3, globals.numSteps);
    rotVel = zeros(3, globals.numSteps);
    rotAcc = zeros(3, globals.numSteps);

    q=zeros(4, globals.numSteps);
    q(:,1)=[0.000001;0.000001;0.000001; sqrt(1-3*.000001^2)];
    qdotSave=zeros(4, globals.numSteps);

    
    %%%% For Kalman Filter %%%%    
    covariance = diag([0.0001 0.00001 0.0000001 0.00001 0.00001 0.0001 0.0001 0.0001 0.0001 0.0001]);    
    % state is a 10x1 matrix of [position velocity quaternions]
    state = [transPos(:,1)' transVel(:,1)' q(:,1)']';
    kalmanHistory = zeros(10,globals.numSteps/kalmanFreq);
    kalmanHistory(:,1)=state;
    kalmanVarHistory=zeros(10,globals.numSteps/kalmanFreq);
    kalmanVarHistory(:,1)=diag(covariance);
    
    
    eBrakesActuated = false;    
    endSimul=[0];
    disp('Simulation Initialized')
    %%%BEGIN LOOPING THROUGH TIMESTEPS%%%
    for n = 2:globals.numSteps
        if sum(endSimul)<10
%         if mod(n,1/(10*globals.timestep)) == 0
%             disp('--------------------------')
%             disp(n*globals.timestep)
%             disp(transPos(:,n-1)')
%             disp(transAcc(:,n-1)')
%         end
        
        %calculate rotation matrix
        q0=q(4,n-1);
        q1=q(1,n-1);
        q2=q(2,n-1);
        q3=q(3,n-1);
        rotMatrix=[1-2*q2^2-2*q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2);...
                   2*(q1*q2+q0*q3) 1-2*q1^2-2*q3^2 2*(q2*q3-q0*q1);...
                   2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) 1-2*q1^2-2*q2^2];

               
            %get all forces in local
            localForces=zeros(3,25);
            localPoints=zeros(3,25);
            forceIndex = 1;

               %%%%% AIR SKATES %%%%%

               % for each air skate, calculate the force
               for i = 1:length(pod.airskate(:,1,1))
                   for j = 1:length(pod.airskate(i,:,:))
                      point= rotMatrix*squeeze(pod.airskate(i,j,:)) + transPos(:,n-1); 
                      vertDist = tube.railHeight + point(3);
                      pointForce=SkateForce(vertDist,50e3,pod.skateSegmentLength);
                        localForces(3,forceIndex)= pointForce/length(pod.airskate(i,:,:));
                        localPoints(:,forceIndex) = squeeze(pod.airskate(i,j,:));
                        forceIndex = forceIndex + 1;

                   end
               end


                %%%%% SPACEX PUSHER %%%%%
                if transPos(1,n-1) < globals.pusherDistance
                    localPusherForce = rotMatrix\[globals.pusherForce; 0; 0];
                    localPusherPoint = pod.COM;

                    localForces(:,forceIndex) = localPusherForce;
                    localPoints(:,forceIndex) = localPusherPoint;
                end
                    forceIndex = forceIndex + 1;
                
                
                %%%%% DAMPING FORCE %%%%%
                localDampForce=[0,0,-transVel(3,n-1)*0.0001];
                localDampPoint=[0;0;0];

                localForces(:,forceIndex) = localDampForce;
                localPoints(:,forceIndex) = localDampPoint;
                forceIndex = forceIndex + 1;
               
                    
                %%%%% LAT CONTROL %%%%%
                localLCFrontForce=-0.05*(dot((rotMatrix*[0.8*pod.length/2;0;0]+transVel(:,n-1)),[0;1;0]))*(abs(dot((rotMatrix*[0.8*pod.length/2;0;0]+transVel(:,n-1)),[0;1;0]))>0.05);
                localLCFrontPoint=[0.8*pod.length/2;0;0];

                localForces(:,forceIndex) = localLCFrontForce;
                localPoints(:,forceIndex) = localLCFrontPoint;
                forceIndex = forceIndex + 1;
                
                localLCBackForce=-0.05*(dot((rotMatrix*[-0.8*pod.length/2;0;0]+transVel(:,n-1)),[0;1;0]))*(abs(dot((rotMatrix*[-0.8*pod.length/2;0;0]+transVel(:,n-1)),[0;1;0]))>0.05);
                localLCBackPoint=[-0.8*pod.length/2;0;0];

                localForces(:,forceIndex) = localLCBackForce;
                localPoints(:,forceIndex) = localLCBackPoint;
                forceIndex = forceIndex + 1;
               

                %%%%% DRAG FORCE %%%%%
                drag = pod.dragCoef*globals.airDensity*pod.height*pod.width/2*(transVel(1,n-1))^2;
                localDragForce = rotMatrix\[-drag;0;0];
                localDragPoint= [pod.length/2; 0; 0];

                localForces(:,forceIndex) = localDragForce;
                localPoints(:,forceIndex) = localDragPoint;
                forceIndex = forceIndex + 1;

                %%%GRAVITY FORCE%%
                gravityForce=[0 0 -1*globals.gravity]* pod.mass;
                localGravityForce=rotMatrix\gravityForce';
                localGravityPoint=pod.COM;

                localForces(:,forceIndex) = localGravityForce;
                localPoints(:,forceIndex) = localGravityPoint;
%                 forceIndex = forceIndex + 1;


            if globals.randomForceNoise

%                 forceSize=size(localForces)
                noise=zeros(3,forceIndex);
                for i=1:forceIndex
                   magForce=norm(localForces(:,i));
                   noise(:,i)=sqrt(globals.noiseModifierConst(i)+globals.noiseModifierLin(i)*abs(magForce-globals.noiseModifierZero(i)))*randn(3,1);%(-1*globals.noiseModifier+2*globals.noiseModifier*rand())*magForce;
                end
                localForces=localForces+noise; 
            end
        
            

        
        %calculate torques in local
        netTorque=[0 0 0];
           for i=1:length(localPoints)
              torque=cross(localPoints(:,i)-pod.COM,localForces(:,i));
              netTorque=netTorque+torque';
           end
                   
%         display(netTorque)

        if (n <= 100)
           netTorque = [0 0 0];            
        end
        
        %get theta accel by tensor\torque
        rotAcc(:,n) = pod.tensor \ transpose(netTorque);
        
        %update omega, calculate omega matrix, update quaternions
        rotVel(:,n) = rotVel(:,n-1) + globals.timestep * rotAcc(:,n);

        OmegaMatrix=[0 rotVel(3,n) -rotVel(2,n) rotVel(1,n);...
                    -rotVel(3,n) 0 rotVel(1,n) rotVel(2,n);...
                     rotVel(2,n) -rotVel(1,n) 0 rotVel(3,n);...
                    -rotVel(1,n) -rotVel(2,n) -rotVel(3,n) 0;];

        qdot = 0.5*OmegaMatrix*q(:,n-1);
        qdotSave(:,n)=qdot;
                
        qn=q(:,n-1)+globals.timestep*qdot;
        normQuat=sqrt(sum(qn.^2));
        qn=qn./normQuat;
        
        
        q(:,n) = qn;

        %calculate roll, pitch, and yaw
        rotPos(:,n-1) = [atan2(2*(q0*q1+q2*q3),1-2*(q1^2 + q2^2)); ...
                         asin(2*(q0*q2 - q3*q1));...
                         atan2(2*(q0*q3+q1*q2),1-2*(q2^2 + q3^2))];

        %get forces in global
        globalForces=rotMatrix*localForces;
        
        netForce=[0 0 0];
        for force=globalForces
            netForce = netForce+force';
        end
        
        if eBrakesActuated
            netForce = netForce + [-pod.eBrakeForce; 0; 0];
        end
        
        if (n < 100)
           netForce = [0 0 0];            
        end
        
        %get accel, velocity, position
        transAcc(:,n)=netForce/pod.mass;
        transVel(:,n) = transVel(:,n-1) + transAcc(:,n)*globals.timestep;
        transPos(:,n) = transPos(:,n-1) + transVel(:,n)*globals.timestep;
        
        %%%%%%%%%%% KALMAN FILTER STEP %%%%%%%%%%%%%%%%%%%%
        
        if mod(n,kalmanFreq) == 0
            
            roll = atan2(2*(q0*q1+ q2*q3),1-2*(q1^2 + q2^2));
            pitch = asin(2*(q0*q2 - q3*q1));
            rolls(n/kalmanFreq) = roll;
            % distance sensors
            distDownData=nan(6,1);
            distRailData=nan(5,1);
            distSideData=nan(5,1);
            distance1Pos = transPos(:,n)*ones(1,6)+rotMatrix*pod.bottomDistancePositions;
            for i=1:6
                [~, distDownData(i)] = DistanceFinder(distance1Pos(:,i)); 
            end

            distance2Pos = transPos(:,n)*ones(1,5)+rotMatrix*pod.downRailDistancePositions;
            for i=1:5
                [distRailData(i)] = TrackDistanceFinder(distance2Pos(:,i)); 
            end

            distance3Pos = transPos(:,n)*ones(1,5)+rotMatrix*pod.sideDistancePositions;
            for i=1:5
                [~, distSideData(i)] = InsideDistance(distance3Pos(:,i)); 
            end
            
%             pitotRead=0.5*globals.airDensity*(dot(transVel(:,n),rotMatrix*pod.pitotDirection))^2;
            photoelectricData=nan(6,3);
            for i=1:3
                photoelectricData(i,1) = photoelectricReadingTop(pod.topPhotoElectricPositions(:,i), state, tube.stripDistances(photoelectricCount(i)+1), pod.peToWall(i), tube,pod);
            end
            for i=1:3
                photoelectricData(i,2) = photoelectricReadingLR(pod.leftPhotoElectricPositions(:,i), state, tube.stripDistances(photoelectricCount(3+i)+1), pod.peToWall(3+i), tube,pod);
            end
            for i=1:3
                photoelectricData(i,3) = photoelectricReadingLR(pod.rightPhotoElectricPositions(:,i), state, tube.stripDistances(photoelectricCount(6+i)+1), pod.peToWall(6+i), tube,pod);
            end
            
            peRecord(:,n/kalmanFreq+1)=[photoelectricData(1:3,1);photoelectricData(1:3,2);photoelectricData(1:3,3)];
            
            [photoelectricCount,photoelectricUse,lastPE,prevPE,crossIn]=photoelectric(photoelectricCount,photoelectricData,lastPE,prevPE,crossIn,n*globals.timestep,pod);
            
            sensorData = [[distDownData]...
                [distRailData; NaN;]...
                [distSideData; NaN;]...
                [nan(6,1)]...%[pitotRead;nan(5,1)]...
                photoelectricData];
            sensorInfo(:,n/kalmanFreq)=[distDownData;distRailData;distSideData;photoelectricData(1:3,1);photoelectricData(1:3,2);photoelectricData(1:3,3);];
            execution =  [0 0 0 0 0 0 0];

%             if mod(n,kalmanFreq*3) == 2
%                 dat1=distDownData'
%                 act1=((rotMatrix(3,:)*pod.bottomDistancePositions(:,:))'+transPos(3,n)+tube.railHeight)
%                 dat2=distRailData'
%                 act2=((rotMatrix(3,:)*pod.downRailDistancePositions(:,:))'+transPos(3,n))'
                
                execution(1) = 1;
                sensorRecord(1:6,n/kalmanFreq)=distDownData;
%             end
%             if mod(n,kalmanFreq*3)==1
                execution(2) = 1;
                sensorRecord(7:11,n/kalmanFreq)=distRailData;
%             end
%             if mod(n,kalmanFreq*100) == 0
%                 cTime=n/10000
%                 cX=transPos(1,n)
%             end
%                 
            
%             if mod(n,kalmanFreq*3) == 0
                execution(3) = 1;
                sensorRecord(12:16,n/kalmanFreq)=distSideData;
%             end
% %             if mod(n,kalmanFreq*20) == 14
% %                 execution(4) = 1;
% %             end
            
%             execution(5:7)=photoelectricUse;
%             sensorRecord(17:25,n/kalmanFreq)=[photoelectricData(1:3,1);ph
%             otoelectricData(1:3,2);photoelectricData(1:3,3);];
            
            IMUData = [transAcc(:,n)' rotVel(:,n)']';
%             disp('imu')
%             disp(IMUData)
            % add random noise
            gyroData(:,n/kalmanFreq)=IMUData(4:6);
  
            if globals.randomIMUNoise
                IMUData = IMUData +...
                    [sqrt(globals.IMUAccelSIMCovConst+globals.IMUAccelSIMCovLin.*abs(IMUData(1:3)-globals.IMUAccelSIMCovZero)).*randn(3,1);...
                    sqrt(globals.IMUGyroSIMCovConst+globals.IMUGyroSIMCovLin.*abs(IMUData(4:6)-globals.IMUGyroSIMCovZero)).*randn(3,1);];
            end
            
            noisyGyroData(:,n/kalmanFreq)=IMUData(4:6);
%           
%             if mod(n,kalmanFreq*5) == 0
%                 dat1=distDownData'
%                 act1=((rotMatrix(3,:)*pod.bottomDistancePositions(:,:))'+transPos(3,n)+tube.railHeight)'
%                 dat2=distRailData'
%                 act2=((rotMatrix(3,:)*pod.downRailDistancePositions(:,:))'+transPos(3,n))'
%             end
            if globals.randomSensorNoise
                sensorData = sensorData + [sqrt(globals.distDownSIMCovConst+globals.distDownSIMCovLin.*abs(distDownData-globals.distDownSIMCovZero)).*randn(6,1)...
                    [sqrt(globals.distDownRailSIMCovConst+globals.distDownRailSIMCovLin.*abs(distRailData-globals.distDownRailSIMCovZero));NaN].*randn(6,1)...
                    [sqrt(globals.distSideSIMCovConst+globals.distSideSIMCovLin.*abs(distSideData-globals.distSideSIMCovZero));NaN].*randn(6,1)...
                    [ nan(6,1)].*randn(6,1)...%[sqrt(globals.pitotSIMCovConst+globals.pitotSIMCovLin.*abs(pitotRead-globals.pitotSIMCovZero)); nan(5,1)].*randn(6,1)...
                    [sqrt(globals.peTopSIMCovConst+globals.peTopSIMCovLin.*abs(photoelectricData(1:3,1)-globals.peTopSIMCovZero)); nan(3,1)]...
                    [sqrt(globals.peLeftSIMCovConst+globals.peLeftSIMCovLin.*abs(photoelectricData(1:3,2)-globals.peLeftSIMCovZero)); nan(3,1)]...
                    [sqrt(globals.peRightSIMCovConst+globals.peRightSIMCovLin.*abs(photoelectricData(1:3,3)-globals.peRightSIMCovZero)); nan(3,1)]];
            end
%              if mod(n,kalmanFreq*5) == 0
%                 dat1=distDownData';
%                 noiseDat1=sensorData(:,1)';
%                 act1=((rotMatrix(3,:)*pod.bottomDistancePositions(:,:))'+transPos(3,n)+tube.railHeight)'
%                   noiseDiff1=act1-noiseDat1
%                   mean(noiseDiff1)
% %                   diff1=act1-dat1
% %                   mean(diff1)
%                 dat2=distRailData';
%                 noiseDat2=sensorData(1:5,2)';
%                 act2=((rotMatrix(3,:)*pod.downRailDistancePositions(:,:))'+transPos(3,n))'
%                   noiseDiff2=act2-noiseDat2
%                   mean(noiseDiff2)
% %                   diff2=act2-dat2
% %                   mean(diff2)
% %  
%             end
            
            
%             if (n/kalmanFreq > 1)
            IMUData = IMUData + [rotMatrix\[0; 0; globals.gravity]; 0;0 ;0;];
            omegaX=IMUData(4);   %angular velocity in x from IMU
            omegaY=IMUData(5);   %angular velocity in y from IMU
            omegaZ=IMUData(6);
            
            q1=noisyKalmanQIMU(1,n/kalmanFreq);  
            q2=noisyKalmanQIMU(2,n/kalmanFreq);  
            q3=noisyKalmanQIMU(3,n/kalmanFreq);  
            q0=noisyKalmanQIMU(4,n/kalmanFreq);  
            noisyKalmanRotIMU=[1-2*q2^2-2*q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2);...
                                 2*(q1*q2+q0*q3) 1-2*q1^2-2*q3^2 2*(q2*q3-q0*q1);...
                                 2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) 1-2*q1^2-2*q2^2;];
            noisyKalmanOmegaIMU=[0 omegaZ -omegaY omegaX;...
                                -omegaZ 0 omegaX omegaY;...
                                omegaY -omegaX 0 omegaZ;...
                                -omegaX -omegaY -omegaZ 0;];
            
            noisyKalmanVelIMU(:,n/kalmanFreq+1) = noisyKalmanVelIMU(:,n/kalmanFreq) + globals.kalmanTimestep *(noisyKalmanRotIMU*(IMUData(1:3))-[0;0;globals.gravity]);
            noisyKalmanPosIMU(:,n/kalmanFreq+1) = noisyKalmanPosIMU(:,n/kalmanFreq) + globals.kalmanTimestep*noisyKalmanVelIMU(:,n/kalmanFreq);
            
            qdotIMU(:,n/kalmanFreq)=0.5*noisyKalmanOmegaIMU*noisyKalmanQIMU(:,n/kalmanFreq);
     

            noisyKalmanQIMU(:,n/kalmanFreq+1)=noisyKalmanQIMU(:,n/kalmanFreq)+0.5*globals.kalmanTimestep*noisyKalmanOmegaIMU*noisyKalmanQIMU(:,n/kalmanFreq);
            
            normQuat=sqrt(sum((noisyKalmanQIMU(:,n/kalmanFreq+1)).^2));
            noisyKalmanQIMU(:,n/kalmanFreq+1)=noisyKalmanQIMU(:,n/kalmanFreq+1)./normQuat;
            %             else
%                 noisyKalmanVelIMU(:,n/kalmanFreq) = globals.kalmanTimestep *(IMUData(1:3));
%             end
            
%                 %disp('vel');
%                 disp(noisyKalmanVelIMU(:,n/kalmanFreq));
%                 %disp('pos');
%                 disp(noisyKalmanPosIMU(:,n/kalmanFreq));                
%                 scannerDistances(n/kalmanFreq) = distDownData;
%             disp(distDownData);
            if globals.sensorFailure
                if rand(1)<globals.failureRate
                    failed=true;
                    st=ceil(7*rand(1));
                    sensorList=sensorData(:,st);
                    if length(sensorList(~isnan(sensorList)))>0
                        sensorType(end+1)=st;
                        sensorNum(end+1)=ceil(rand(1)*length(find(~isnan(sensorList))));
                        if rand(1)>0.5
                            sensorFailState(end+1)=1;
                            sensorData(sensorNum(end),sensorType(end))=globals.sensorMaxs(sensorNum(end),sensorType(end));
                            fprintf('*****Sensor %d of type %d has failed at %d milliseconds, reading high.********\n',sensorNum(end),sensorType(end),n*globals.timestep*1000)
                        else
                            sensorFailState(end+1)=0;
                            sensorData(sensorNum(end),sensorType(end))=globals.sensorMins(sensorNum(end),sensorType(end));
                            fprintf('*****Sensor %d of type %d has failed at %d milliseconds, reading low.********\n',sensorNum(end),sensorType(end),n*globals.timestep*1000)
                        end
                    end
                end
                if globals.failurePersist && failed
                    for i=1:length(sensorType)
                       if sensorFailState(i)==1
                           sensorData(sensorNum(i),sensorType(i))=globals.sensorMaxs(sensorNum(i),sensorType(i));
                       elseif sensorFailState(i)==0
                           sensorData(sensorNum(i),sensorType(i))=globals.sensorMins(sensorNum(i),sensorType(i));
                       else
                           fprintf('**********************SOMETHING WRONG**********************\n')
                       end
                    end
                end
            end

            [sensorUse,numberUsed] = sensorFailureDetection(sensorData,globals,pod,tube);
           
            [state, covariance,~] = KalmanFilterHyperloop(state, covariance, IMUData, sensorData, execution, sensorUse,numberUsed, n/kalmanFreq,photoelectricCount,globals,pod,tube);
%              if mod(n,kalmanFreq*5) == 0
%                 distance1Pos
%                 transPos(:,n)
%                 state
%                 noisyKalmanPosIMU
%                 noisyKalmanPosIMU+pod.height/2
%             end
            
            kalmanHistory(:,n/kalmanFreq+1) = state;
            kalmanVarHistory(:,n/kalmanFreq+1) = diag(covariance);
%             disp('---------------------')
%             disp(n/kalmanFreq)
%             disp(kalmanHistory(3,n/kalmanFreq))
%             disp(noisyKalmanPosIMU(3,n/kalmanFreq))
%             disp(transPos(3,n))
        end
        
        
%         %check collisions
%         endSimul=false;
       for point = pod.collisionPoints
            realPoint = rotMatrix*point + transPos(:,n);
            [distance, ~]= DistanceFinder(realPoint);
            if mod(n,kalmanFreq) == 0
                q0=kalmanHistory(10, n/kalmanFreq);
                q1=kalmanHistory(7, n/kalmanFreq);
                q2=kalmanHistory(8, n/kalmanFreq);
                q3=kalmanHistory(9, n/kalmanFreq);
                kalmanRotMatrix=[1-2*q2^2-2*q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2);...
                           2*(q1*q2+q0*q3) 1-2*q1^2-2*q3^2 2*(q2*q3-q0*q1);...
                           2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) 1-2*q1^2-2*q2^2];
                kalmanPoint = kalmanRotMatrix*point + kalmanHistory(1:3,n/kalmanFreq);
                [kalmanDistance, ~] = DistanceFinder(kalmanPoint);
                if kalmanDistance <= 0
                   disp('simulation Ended due to Kalman thinking we crashed');
                   disp(point);
                   disp(realPoint);
                   disp(kalmanPoint);
                   disp(distance);
                   disp(kalmanDistance);
%                    endSimul(end+1)=1;
                   break;
                end
            end

            if distance<=0
                norm(point)
                disp('Collision Occurred at timestep');
                disp(n);
%                 endSimul(end+1)=1;
                   disp(point);
                   disp(realPoint);
                   disp(kalmanPoint);
                   disp(distance);
                   disp(kalmanDistance);
                break
            end
        end
        end
    end
    
    display(size(kalmanHistory(3,:)))
    display(size(transPos(3,:)))
    figure;
    subplot(3,1,1)
    hold all
    h1=plot(kalmanHistory(1,:),'b');
    h1a=plot(kalmanHistory(1,:)-sqrt(kalmanVarHistory(1,:)),'b:');
    h1b=plot(kalmanHistory(1,:)+sqrt(kalmanVarHistory(1,:)),'b:');
    h2=plot(downsample(transPos(1,:),kalmanFreq),'g');
    h3=plot(noisyKalmanPosIMU(1,:),'r');
    legend([h1 h2 h3],{'Kalman','Physical','IMU'})
    title('X-position')
    hold off
    subplot(3,1,2)
    hold all
    h1=plot(kalmanHistory(2,:),'b');
    h1a=plot(kalmanHistory(2,:)-sqrt(kalmanVarHistory(2,:)),'b:');
    h1b=plot(kalmanHistory(2,:)+sqrt(kalmanVarHistory(2,:)),'b:');
    h2=plot(downsample(transPos(2,:),kalmanFreq),'g');
    h3=plot(noisyKalmanPosIMU(2,:),'r');
    legend([h1 h2 h3],{'Kalman','Physical','IMU'})
    title('Y-position')
    for i=12:16
        plot(sensorInfo(i,:),'k:')
%         plot(sensorRecord(i,:),'k-')
    end
    hold off
    subplot(3,1,3)
    hold all
    h1=plot(kalmanHistory(3,:),'b');
    h1a=plot(kalmanHistory(3,:)-sqrt(kalmanVarHistory(3,:)),'b:');
    h1b=plot(kalmanHistory(3,:)+sqrt(kalmanVarHistory(3,:)),'b:');
    h2=plot(downsample(transPos(3,:),kalmanFreq),'g');
    h3=plot(noisyKalmanPosIMU(3,:),'r');
    legend([h1 h2 h3],{'Kalman','Physical','IMU'})
    for i=1:6
        plot(sensorInfo(i,:),'k:')
%         plot(sensorRecord(i,:),'k-')
    end
    for i=7:11
        plot(sensorInfo(i,:),'c:')
%         plot(sensorRecord(i,:),'c-')
    end
    title('Z-position')
    hold off
    
    
    
    figure;
    subplot(3,1,3)
    hold all
    IMUvReal=noisyKalmanVelIMU(3,1:end-1)-downsample(transVel(3,:),kalmanFreq);
    plot(kalmanHistory(6,:),'b');
    plot(kalmanHistory(6,:)-sqrt(kalmanVarHistory(6,:)),'b:');
    plot(kalmanHistory(6,:)+sqrt(kalmanVarHistory(6,:)),'b:');
    plot(downsample(transVel(3,:),kalmanFreq),'g');
    plot(noisyKalmanVelIMU(3,:),'r');
    title('Z-velocity')
    
    subplot(3,1,1)
    hold all
%     IMUvReal=noisyKalmanVelIMU(3,1:end-1)-downsample(transVel(3,:),kalman
%     Freq);
    plot(kalmanHistory(4,:),'b')
    plot(kalmanHistory(4,:)-sqrt(kalmanVarHistory(4,:)),'b:');
    plot(kalmanHistory(4,:)+sqrt(kalmanVarHistory(4,:)),'b:');
    plot(downsample(transVel(1,:),kalmanFreq),'g');
    plot(noisyKalmanVelIMU(1,:),'r');
    title('X-velocity')
    subplot(3,1,2)
    hold all
%     IMUvReal=noisyKalmanVelIMU(3,1:end-1)-downsample(transVel(3,:),kalmanFreq);

    h1=plot(kalmanHistory(5,:),'b');
    h1a=plot(kalmanHistory(5,:)-sqrt(kalmanVarHistory(5,:)),'b:');
    h1b=plot(kalmanHistory(5,:)+sqrt(kalmanVarHistory(5,:)),'b:');
    h2=plot(downsample(transVel(2,:),kalmanFreq),'g');
    h3=plot(noisyKalmanVelIMU(2,:),'r');
    title('Y-velocity')
    
    legend([h1 h2 h3],{'Kalman','Physical','IMU'})
    hold off
    
    
    figure;
    subplot(4,1,1)
    hold all
    plot(kalmanHistory(7,:),'b');
    plot(kalmanHistory(7,:)-sqrt(kalmanVarHistory(7,:)),'b:');
    plot(kalmanHistory(7,:)+sqrt(kalmanVarHistory(7,:)),'b:');
    plot(downsample(q(1,:),kalmanFreq),'g');
    plot(noisyKalmanQIMU(1,:),'r');
    title('q1')
    subplot(4,1,2)
    hold all
%     imuvreal=noisykalmanvelimu(3,1:end-1)-downsample(transvel(3,:),kalman
%     freq);
    plot(kalmanHistory(8,:),'b');
    plot(kalmanHistory(8,:)-sqrt(kalmanVarHistory(8,:)),'b:');
    plot(kalmanHistory(8,:)+sqrt(kalmanVarHistory(8,:)),'b:');
    plot(downsample(q(2,:),kalmanFreq),'g');
    plot(noisyKalmanQIMU(2,:),'r');
    title('q2')
    subplot(4,1,3)
    hold all
%     imuvreal=noisykalmanvelimu(3,1:end-1)-downsample(transvel(3,:),kalmanFreq);

  plot(kalmanHistory(9,:),'b');
    plot(kalmanHistory(9,:)-sqrt(kalmanVarHistory(9,:)),'b:');
    plot(kalmanHistory(9,:)+sqrt(kalmanVarHistory(9,:)),'b:');
    plot(downsample(q(3,:),kalmanFreq),'g');
    plot(noisyKalmanQIMU(3,:),'r');
    title('q3')
    subplot(4,1,4)
    hold all
    h1=plot(kalmanHistory(10,:),'b');
    h1a=plot(kalmanHistory(10,:)-sqrt(kalmanVarHistory(10,:)),'b:');
    h1b=plot(kalmanHistory(10,:)+sqrt(kalmanVarHistory(10,:)),'b:');
    h2=plot(downsample(q(4,:),kalmanFreq),'g');
    h3=plot(noisyKalmanQIMU(4,:),'r');
    title('q0')
    
    legend([h1 h2 h3],{'kalman','physical','imu'})
    hold off

    kalmanRoll = atan2(2*(kalmanHistory(10,:).*kalmanHistory(7,:)+ kalmanHistory(8,:).*kalmanHistory(9,:)),1-2*(kalmanHistory(7,:).^2 + kalmanHistory(8,:).^2));
            kalmanPitch = asin(2*(kalmanHistory(10,:).*kalmanHistory(8,:) - kalmanHistory(9,:).*kalmanHistory(7,:)));
            kalmanYaw= atan2(2*(kalmanHistory(10,:).*kalmanHistory(9,:)+kalmanHistory(7,:).*kalmanHistory(8,:)),1-2*(kalmanHistory(8,:).^2 + kalmanHistory(9,:).^2));
    
            
            imuRoll = atan2(2*(noisyKalmanQIMU(4,:).*noisyKalmanQIMU(1,:)+ noisyKalmanQIMU(2,:).*noisyKalmanQIMU(3,:)),1-2*(noisyKalmanQIMU(1,:).^2 + noisyKalmanQIMU(2,:).^2));
            imuPitch = asin(2*(noisyKalmanQIMU(4,:).*noisyKalmanQIMU(2,:) - noisyKalmanQIMU(3,:).*noisyKalmanQIMU(1,:)));
            imuYaw= atan2(2*(noisyKalmanQIMU(4,:).*noisyKalmanQIMU(3,:)+noisyKalmanQIMU(1,:).*noisyKalmanQIMU(2,:)),1-2*(noisyKalmanQIMU(2,:).^2 + noisyKalmanQIMU(3,:).^2));
    
           
            
figure;
    subplot(3,1,3)
    hold all
    plot(kalmanRoll,'b');
    plot(downsample(rotPos(1,:),kalmanFreq),'g');
    plot(imuRoll,'r');
    title('roll')
    subplot(3,1,1)
    hold all
    plot(kalmanPitch,'b')
    plot(downsample(rotPos(2,:),kalmanFreq),'g');
    plot(imuPitch,'r');
    title('pitch')
    subplot(3,1,2)
    hold all
    h1=plot(kalmanYaw,'b');
    h2=plot(downsample(rotPos(3,:),kalmanFreq),'g');
    h3=plot(imuYaw,'r');
    title('yaw')
    
    legend([h1 h2 h3],{'kalman','physical','imu'})
    hold off
   


% end
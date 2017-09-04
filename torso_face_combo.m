%% Copyright 2012 The MathWorks, Inc.
% Create a detector object

srcFiles = dir('C:\Users\SAURADIP\Documents\MATLAB\video frame multimodal\new race fused 3\seed\*.jpg');  % the folder in which ur images exists
for g = 22 : 22
  
    filename = strcat('C:\Users\SAURADIP\Documents\MATLAB\video frame multimodal\new race fused 3\seed\',srcFiles(g).name);
   
faceDetector = vision.CascadeObjectDetector;   

% Read input image
I = imread(filename);

% Detect faces
bbox = step(faceDetector, I); 

% Create a shape inserter object to draw bounding boxes around detections
shapeInserter = vision.ShapeInserter('BorderColor','Custom','CustomBorderColor',[255 255 0]); 

% Draw boxes around detected faces and display results              
I_faces = step(shapeInserter, I, int32(bbox));    

%%
% Create a detector object 
faceDetector = vision.CascadeObjectDetector; 
faceDetector.MinSize = [40 40];
faceDetector.MaxSize = [80 80];  % fixed condition for a marathon athlete too close to camera , face detected for minsize [30 30] and maxsize[80 80] and scaleFactor 1.014 
faceDetector.ScaleFactor = 1.015;


bodyDetector = vision.CascadeObjectDetector('UpperBody'); 
bodyDetector.MinSize = [220 220]; 
bodyDetector.MaxSize = [280 280];% fixed condition for a marathon athlete too close to camera , torso detected for minsize [130 130] and scaleFactor 1.01 
bodyDetector.ScaleFactor = 1.0051;

% Read input image and detect upper body
bbox_body = step(bodyDetector, I);
% Read input image and detect upper body
bbox_face = step(faceDetector, I);

% Draw bounding boxes
shapeInserter = vision.ShapeInserter('BorderColor','Custom','CustomBorderColor',[255 255 0]);
I_body = step(shapeInserter, I, int32(bbox_face));
figure, imshow(I_body);
%%
% bbox_face = zeros(size(bbox_body));
% if(size(bbox_face,1)>0)
% for i=1:size(bbox_face,1)
%     Icrop = imcrop(I,bbox_face(i,:));
%     bbox = step(faceDetector,Icrop);
%     %bbox_face(i,:) = bbox + [bbox_body(i,1:2)-1 0 0];
% end
% end
%     
% if(size(bbox_body,1)>0)
% for i=1:size(bbox_body,1)
%     Icrop2 = imcrop(I,bbox_body(i,:));
%     bbox = step(bodyDetector,Icrop2);
%     %bbox_face(i,:) = bbox + [bbox_body(i,1:2)-1 0 0];
% end
% end
% I_faces2 = step(shapeInserter, I, int32(bbox_face));
% figure, imshow(I_faces2);

img = I;
% 2) Display the image and store its handle:
h_im = imshow(img);
% 3) Create an ellipse defining a ROI:
temp = cat(3, zeros(size(I,1),size(I,2)),zeros(size(I,1),size(I,2)),zeros(size(I,1),size(I,2)));
temp1=temp;


%figure, imshow(ROI);

folder = strcat('frame', num2str(g)); 

[status,msg]=mkdir('C:\Users\SAURADIP\Documents\MATLAB\video frame multimodal\new race fused 3\detection result\',folder);
path=strcat('C:\Users\SAURADIP\Documents\MATLAB\video frame multimodal\new race fused 3\detection result\',folder);
for i=1:size(bbox_face,1)
e = imrect(gca,bbox_face(i,:));
%f= imrect(gca,bbox_body(i,:));
% 4) Create a mask from the ellipse:
BW = createMask(e,h_im);

% 4a) (For color images only) Augment the mask to three channels:
BW(:,:,2) = BW;
BW(:,:,3) = BW(:,:,1);

temp=temp + BW;
ROI = img;  % face only 
ROI(temp == 0) = 0;
% torso = fullfile(path, sprintf('combo_image_%05d.JPG', i));   %Writing the final seed frame
%                  imwrite(im2bw(ROI), torso);
% torsocolor = fullfile(path, sprintf('combo_image_colored%05d.JPG', i));   %Writing the final seed frame
%                  imwrite(ROI, torsocolor);
end

for i=1:size(bbox_body,1)
    
f= imrect(gca,bbox_body(i,:));
BW2 = createMask(f,h_im);
BW2(:,:,2) = BW2;
BW2(:,:,3) = BW2(:,:,1);
temp1=temp1+BW2;
ROI1 = img; % body only
ROI1(temp1 == 0) = 0;
end

temp2=ROI1+ROI;  % face + body only


ROI3 = img;
ROI3(temp2 == 0) = 0;
torsofinal = fullfile(path, sprintf('combo_binary%05d.JPG',''));   %Writing the final seed frame
                 imwrite(im2bw(ROI3), torsofinal);
                 torsofinal = fullfile(path, sprintf('combo_color%05d.JPG','')); 
                 imwrite(ROI3,torsofinal);


end


%         points = detectMinEigenFeatures(rgb2gray(ROI3),'ROI',bbox_face(1,:));
% 
%         pointImage = insertMarker(ROI3,points.Location,'+','Color','white');
% figure;
% imshow(pointImage);
%          tracker = vision.PointTracker('MaxBidirectionalError',1);
%           pointTracker = vision.PointTracker;
%          initialize(tracker,points.Location,I);
%           [new_points,point_validity,scores] = step(tracker,I);
%         
%         
%   

%% Copyright 2012 The MathWorks, Inc.
% Create a detector object

srcFiles = dir('C:\Users\SAURADIP\Documents\MATLAB\video frame multimodal\new race fused 3\seed\*.jpg');  % the folder in which ur images exists
for g = 9 : 9
  
    filename = strcat('C:\Users\SAURADIP\Documents\MATLAB\video frame multimodal\new race fused 3\seed\',srcFiles(g).name);
    %color = strcat('C:\Users\SAURADIP\Documents\MATLAB\video frame multimodal\new race fused 3\color\',srcFiles(g).name);
   
faceDetector = vision.CascadeObjectDetector;
% colorskin = imread(color);
%  cropped=imcrop(colorskin,[ 0 0 141 141]);
%    resize=imresize(cropped,4);
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
faceDetector.ScaleFactor = 1.015;  % original value is 1.015 changed to 1.016


bodyDetector = vision.CascadeObjectDetector('UpperBody'); 
bodyDetector.MinSize = [160 160]; % [ 160 160 ] for case 4 : no face but torso
bodyDetector.MaxSize = [280 280];% fixed condition for a marathon athlete too close to camera , torso detected for minsize [130 130] and scaleFactor 1.01 
bodyDetector.ScaleFactor = 1.006; %  for case 4 : no face but torso 1.006 case gives good result 

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
ROI=0;
ROI1=0;
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
bod=ROI1;
bod1=ROI1(temp1 == 0);
temp2=ROI1+ROI; % face + body only
ROI3 = img;
ROI3(temp2 == 0) = 0;

temp3=im2bw(temp+temp1); % face + body only but binary rectangle form

stats = regionprops(temp3,'BoundingBox'); % gives the rectanglestarting pixel width height from binary image of face + body , contains list of all rectangles/regions in image
valid = zeros(10,1);
changed=zeros(10);
binary = im2bw(ROI3);
flag=0;
face=zeros(10,4);
body=zeros(10,4);
for m=1:size(stats,1)
    M = stats(m).BoundingBox;
    if ( (M(1,3) <= 120) && (M(1,4) <= 120 ))
        
        face(m,:)=M;
    else
        
        body(m,:)=M;
    end
end
face( all(~face,2), : ) = [];
body( all(~body,2), : ) = [];
        for k=1:size(face,1)
           overlapRatio=0; % for case 1 
           
               
%               temp_face= face(k,4);   % for case 3
%               face(k,3)= 2*(face(k,3)); % for case 3
%               face(k,1)= (face(k,1)-20); % for case 3
              
              
%             box=bbox_face(k,:);
%             x_axis=box(1,1);
%             y_axis=box(1,2);
%             width=box(1,3);
%             height=box(1,4);
           
           if(face(k,1) ~=0 || face(k,2) ~=0)
                         limit=0; 
                        %overlapRatio=0; % for case 1
                     for j=1:size(body,1)  
                         temp_body=body(k,4); % case 4 
                         h = temp_body / 4; % case 4 
                         tot_limit= 7*h ; % case 4 
                         limit = tot_limit-(body(j,4)); % case 4
                        % limit = 564-(face(k,4)+face(k,2)); % for case 3
                         
                         % for case 1 below while()
                         %while ( overlapRatio == 0  &&  564-(face(k,4)+face(k,2))>=0 )
                         %while ( limit >=0 && (face(k,4)) <= 7*temp_face) % for case 3
                         while ( limit >=0 ) % case 4
                        changed(k,1)=face(k,4); % for all
                            face(k,4)=face(k,4)+1; % for case 3 and 1
                            body(j,4)=body(j,4)-1;
                            body(j,2)=body(j,2)-1;% for case 4 
                            body(j,4)= body(j,4)+1; % case 4
                            limit = tot_limit-(body(j,4)); %for case 4
                            %limit = 564-(face(k,4)+face(k,2)); % for case 3
                             overlapRatio =rectint(face(k,:),body(j,:)); % for case 1
                           
                        
                        end
                     end
                
                 
                   
           end
                    
                
            
             
        end
        
        
        
   % for all cases      
        
ROI=0;    
temp=0;
BW=0;
ROI1=0;
for i=1:size(face,1)
e = imrect(gca,face(i,:));
%f= imrect(gca,bbox_body(i,:));
% 4) Create a mask from the ellipse:
BW = createMask(e,h_im);

% 4a) (For color images only) Augment the mask to three channels:
BW(:,:,2) = BW;
BW(:,:,3) = BW(:,:,1);

temp=temp + BW;
ROI = img;  % face only 
ROI(temp == 0) = 0;   
end

for i=1:size(body,1)
    
f= imrect(gca,body(i,:));
BW2 = createMask(f,h_im);
BW2(:,:,2) = BW2;
BW2(:,:,3) = BW2(:,:,1);
temp1=temp1+BW2;
ROI1 = img; % body only
ROI1(temp1 == 0) = 0;
end

temp2=ROI1+ROI; % face + body only
ROI3 = img;
temp4=temp1+temp;
ROI3(temp2 == 0) = 0;

torsofinal = fullfile(path, sprintf('connect_binary%05d.JPG',''));   %Writing the final seed frame
                 imwrite(im2bw(ROI3), torsofinal);
                 torsofinal = fullfile(path, sprintf('connect_color%05d.JPG','')); 
                 imwrite(ROI3,torsofinal);

   
                 
                 
 % for text detection we need to find k=2 on fused image where input image
 % is the corrected human body with face
        
%resize=imresize(bod,0.5);           
% B=rgb2gray(resize);
%                 
% C=double(B)/double(max(max(B)));
%  [cr,cc]=size(C);
% %DC=zeros(cr+3,cc+3);
% D=double(C);
% 
% DC=D;
% win=5;  % sobel window size , generaly 3 
% 
%  
% switch win
%     case 3
%        Sx=[-1 -2 -1 ; 0 0 0 ; 1 2 1]; %for 3 x 3
%        Sy=[-1 0 1 ; -2 0 2 ; -1 0 1];
%     case 7 
%     Sx = [ 3 2 1   0 -1 -2 -3
%            4 3  2  0 -2 -3 -4 
%            5 4  3  0 -3 -4 -5 
%            6 5  4  0 -4 -5 -6 
%            5 4  3  0 -3 -4 -5 
%            4 3  2  0 -2 -3 -4
%            3 2  1  0 -1 -2 -3 ];
%        Sy=transpose(Sx);
%     case 5
%        Sx=[-1  -4  -6  -4  -1 ; -2  -8  -12  -8  -2 ; 0  0  0  0  0  ; 2  8  12  8  2 ; 1  4  6  4  1 ];
%        Sy=[1  2   0   -2  -1 ; 4   8   0   -8  -4 ; 6   12   0   -12  -6 ; 4   8   0   -8  -4 ; 1   2   0   -2  -1];
% end 
%         [r,c]=size(D);
%         
%         O=zeros(r,c);
%         O1=zeros(r,c);
%        
%         f=zeros(r-win,c-win);
%         m=zeros(r-win,c-win);
%          p=zeros(r-win,c-win);
%         px=zeros(win,win);
%          py=zeros(win,win);
%          xyz=zeros(r-win,c-win);
%          xy=zeros(r-win,c-win);
%           pic=zeros(r-win,c-win);
%           fuse=zeros(r-win,c-win);
%           %MG_dir=zeros(r-win,c-win);
%           DC_dir=zeros(r-win,c-win);
%           dc=xyz;
%           dc_lam1=zeros(r-win,c-win);
%           dc_lam2=zeros(r-win,c-win);
%           mg_division=DC_dir;
%           dc_mag=DC_dir;
%           DC_quad=DC_dir;
%           MG_quad=MG_dir;
%           mag=m;
%           %means=mag; % for which i need to find k means cluster
%           if(r>0 )
%             for k=1:(r-win)
%             for l=1:(c-win)
%         patch=D(k:(k+(win-1)),l:(l+(win-1)));
%         gx=Sx.*patch;
%         gy=Sy.*patch;
%        Gx=sum(sum(gx));
%        Gy=sum(sum(gy));
%         O(k,l)=sum(sum(gx));
%         O1(k,l)=sum(sum(gy));
% %             for i=k:(k+2)
% %             for j=l:(l+2)
% %       O(i,j)=sum(sum(gx));
%        % O1(i,j)=sum(sum(gy));
% %             end
% %             end
%         sums=0.0;
%        
%        
%   
%        %m(k,l)=log(sum(f(:))/double(9));
%         
%             end
%             end
%  f=power(sqrt((power(O,2)+power(O1,2))/2),4); 
%  Ix_binary=round(double(O)/double(max(max(O))),1);
%  Iy_binary=round(double(O1)/double(max(max(O1))),1);
%    for kk=1:(r-win)
%        for ll=1:(c-win)
%         G=f(kk:(kk+(win-1)),ll:(ll+(win-1)));
%         if ( (sum(sum(G))/double(win*win))~=0)
%             if(log(sum(sum(G))/double(win*win)) > 0 )
%                m(kk,ll)=((log1p(sum(sum(G))/double(win*win))));   % Magnitude Gradient of Image
%             else
%             m(kk,ll)= abs((log1p(sum(sum(G))/double(win*win))));
%             end
%         end
%         if(Iy_binary(kk,ll)~=0)
%         mg_division(kk,ll)=round((Ix_binary(kk,ll)./Iy_binary(kk,ll)),2); % division of Iy / Ix for each pixel of image
%         
%         end
% px=(O(kk:(kk+(win-1)),ll:(ll+(win-1))));
% py=(O1(kk:(kk+(win-1)),ll:(ll+(win-1))));
%       T=cov(px,py);
%       Ei=eig(T);
%       dc_lam1(kk,ll)=Ei(1,1);
%       dc_lam2(kk,ll)=Ei(2,1);
%       
%        x1=Ei(1,1);
%      y1=Ei(2,1);
%     %DC_dir(kk,ll) =atan(x1/y1); % direction of each edge point for DC
% %         if(x1<0)
% %          x1=0;
% %         end
% %         if(y1<0)
% %          y1=0;
% %         end
%      su=abs(x1-y1);
%      ad=abs(x1+y1);
%      if( dc_lam2(kk,ll)~=0)
%          dc_mag(kk,ll)=round((dc_lam1(kk,ll)./dc_lam2(kk,ll)),1); 
%      end
%      if(ad>0 && su < ad )
%      dc(kk,ll)=power(round(su/ad,2),2);
%      end
%      %xyz(kk,ll)=round(su.*su,4);
%      %xy(kk,ll)=round(ad.*ad,4);
%       if(ad > 0  && (dc(kk,ll)>0))
%      pic(kk,ll)=abs(dc(kk,ll));  % Directional Coherence for the image
%       else
%           pic(kk,ll)=0;
%        
%       end
%       
%       
%        end
%    end
%    
%    mg_division_normalised = round(double(mg_division)/double(max(max(mg_division))),4);  % division of Iy / Ix for each pixel of image and normalised within 0 to 1 
%         MG_dir =atand(mg_division);  % direction of each edge point for MG by tan-1(Iy/Ix)
%    
%    %dc_mag=round(abs(dc_lam1./dc_lam2),1);
%    dc_mag_normalised = round(double(dc_mag)/double(max(max(dc_mag))),1);  % normalised within 0 and 1 result of lambda1/lambda2 rounded to 1 decimal place 
%    DC_dir = atand(dc_mag);  % direction of DC values in DC matrix tan-1(lambda1/lambda2)
%    
%    mag=double(m)/double(max(max(m)));        
%      %dc=double(pic)/double(max(max(pic)));
%      [rr,cc]=size(m);
%      
%      % FOR 5 X 5 Window
% %          for is=1:(rr-5)
% %             for js=1:(cc-5)
% %         matmg = m(is:(is+4),js:(js+4));
% %         matdc = pic(is:(is+4),js:(js+4));
% %         maxmg= max(matmg(:));
% %         minmg= min(matmg(:));
% %          maxdc= max(matdc(:));
% %         mindc= min(matdc(:));
% %         diffmg=maxmg - minmg;
% %         diffdc=maxdc - mindc;
% % %         if(diffmg > diffdc)
% % %             result=diffmg;
% % %         else
% % %             result=diffdc;
% % %         end
% %         %fuse(is,js)=result;
% % %                         if ((diffmg+diffdc)<0)
% % %              fuse(is,js)=0;
% % %          
% % %                         else
% %              fuse(is,js)=diffmg+diffdc;
% %              newfuse=double(fuse)/double(max(max(fuse)));
% %                       %  end            
% %            end
% %          end
%          mg1=imresize(imgradient(im2bw(bod)),.25);
%       fuse=(mag+pic);
%       
%       %fuse=double(fuse)/double(max(max(fuse)));         
%           end
%            means=pic; % for which i need to find k means cluster
% nrows = size(means,1);
% ncols = size(means,2);
% ab = reshape(means,nrows*ncols,1);
% [idx3,centroid] = kmeans(ab,2);
% 
% 
% % Average / Max  / Min cluster will depend on the centroid value , higher
% % centroid means max , lower centroid means min and rest is avg
% 
% cluster1=ab(idx3==1);
% cluster2=ab(idx3==2);
% %cluster3=ab(idx3==3);
% [r1,c1]=size(cluster1); 
% [r2,c2]=size(cluster2); 
% %[r3,c3]=size(cluster3);
% 
% 
% 
% [rf,cf]=size(means);
% df3=zeros(rf,cf);
% df2=zeros(rf,cf);
% df1=zeros(rf,cf);
% df4=zeros(rf,cf);
% 
% %converting the clustered array into real clustered matrix of m x n size
%  
%  for ccol=1:r1
%      % process 1 %
%         value=round(cluster1(ccol,1),2);
%        df1(round(means,2)==value)=1;        %max cluster
%  end
%       
%  for ccol=1:r2
%        value1=round(cluster2(ccol,1),2);
%        df2(round(means,2)==value1)=1;       %avg cluster
%        
%  end
%  
% %  for ccol=1:r3
% %        value2=round(cluster3(ccol,1),4);
% %        df3(round(fuse,4)==value2)=1;        %min cluster
% %       
% %         % process 2 %
% %        
% %  end
%  
%  for i=1:rf
%      for j=1:cf
%        if (df1(i,j) ~=1)
%            df1(i,j)=0;
%        end
%        if (df2(i,j) ~=1)
%            df2(i,j)=0;
%        end
% %         if(df2(i,j))
% %             df4(i,j)=D(i,j);
% %         end
% %        if (df3(i,j) ~=1)
% %            df3(i,j)=0;
% %        end
%            
%      end
%  end
%                 
%                  
%                  
%                  
%       df3=imresize(df1,2);
%  df4=imresize(df2,2);            
%                  
%  x_axis_t=190;
%  y_axis_t=102;
%  
%  x_axis_nt=138;
%  y_axis_nt=102;
%  
%  x_axis_bg=79;
%  y_axis_bg=48;
%  
%  % Direction of Magnitude Gradient 
%  
%  text_MG = 180.*MG_dir(x_axis_t:x_axis_t+7,y_axis_t:y_axis_t+7);
%  nontext_MG = 180.*MG_dir(x_axis_nt:x_axis_nt+7,y_axis_nt:y_axis_nt+7);
%  background_MG = 180.*MG_dir(x_axis_bg:x_axis_bg+7,y_axis_bg:y_axis_bg+7);
%  
%  % Direction of Directional Coherence
%  
%  text_DC = 180.*DC_dir(x_axis_t:x_axis_t+7,y_axis_t:y_axis_t+7);
%  nontext_DC = 180.*DC_dir(x_axis_nt:x_axis_nt+7,y_axis_nt:y_axis_nt+7);
%  background_DC = 180.*DC_dir(x_axis_bg:x_axis_bg+7,y_axis_bg:y_axis_bg+7);
% 
%  
%  % Normal Magnitude Gradient value
%  
%  text_MG_normal = mag(x_axis_t:x_axis_t+7,y_axis_t:y_axis_t+7);
%  nontext_MG_normal = mag(x_axis_nt:x_axis_nt+7,y_axis_nt:y_axis_nt+7);
%  background_MG_normal = mag(x_axis_bg:x_axis_bg+7,y_axis_bg:y_axis_bg+7);
%  
%  % Normal DC Value
%  
%   text_DC_normal = pic(x_axis_t:x_axis_t+7,y_axis_t:y_axis_t+7);
%  nontext_DC_normal = pic(x_axis_nt:x_axis_nt+7,y_axis_nt:y_axis_nt+7);
%  background_DC_normal = pic(x_axis_bg:x_axis_bg+7,y_axis_bg:y_axis_bg+7);
%  
%  
%  
%  % Pixel information 
%  
%  text_pixel = D(x_axis_t:x_axis_t+7,y_axis_t:y_axis_t+7);
%  nontext_pixel = D(x_axis_nt:x_axis_nt+7,y_axis_nt:y_axis_nt+7);
%  background_pixel = D(x_axis_bg:x_axis_bg+7,y_axis_bg:y_axis_bg+7);
%  
%  
%  for l = 1: kk
%      for m =1 : ll
%          if(MG_dir(l,m)>0)
%              MG_dir(l,m)=round(MG_dir(l,m),2);
%              DC_dir(l,m)=round(DC_dir(l,m),2);
%          else
%              MG_dir(l,m)=round(360-MG_dir(l,m),2);
%              DC_dir(l,m)=round(360-DC_dir(l,m),2);
%          end
%      end
%  end
%  
%  
%  for xa=1:kk
%      for ya=1:ll
%          if((DC_dir(xa,ya)<=45.00 && DC_dir(xa,ya)>=0.00) || ( MG_dir(xa,ya)<=45.00 && MG_dir(xa,ya)>=0.00))
%              DC_quad(xa,ya)=1;
%              MG_quad(xa,ya)=1;
%          elseif((DC_dir(xa,ya)<=90.00 && DC_dir(xa,ya)>45.00 ) || (MG_dir(xa,ya)<=45.00 && MG_dir(xa,ya)>45.00))
%                 DC_quad(xa,ya)=2;
%                 MG_quad(xa,ya)=2;
%          elseif((DC_dir(xa,ya)<=135.00 && DC_dir(xa,ya)>90.00 ) || (MG_dir(xa,ya)<=135.00 && MG_dir(xa,ya)>90.00))  
%                 DC_quad(xa,ya)=3;
%                 MG_quad(xa,ya)=3;
%          elseif((DC_dir(xa,ya)<=180.00 && DC_dir(xa,ya)>135.00 ) || (MG_dir(xa,ya)<=180.00 && MG_dir(xa,ya)>135.00))  
%                 DC_quad(xa,ya)=4;
%                 MG_quad(xa,ya)=4;
%          elseif((DC_dir(xa,ya)<=225.00 && DC_dir(xa,ya)>180.00 ) || (MG_dir(xa,ya)<=225.00 && MG_dir(xa,ya)>180.00))  
%                 DC_quad(xa,ya)=5;
%                 MG_quad(xa,ya)=5;
%          elseif((DC_dir(xa,ya)<=270.00 && DC_dir(xa,ya)>225.00 ) || (MG_dir(xa,ya)<=270.00 && MG_dir(xa,ya)>225.00))  
%                 DC_quad(xa,ya)=6;
%                 MG_quad(xa,ya)=6;
%           elseif((DC_dir(xa,ya)<=315.00 && DC_dir(xa,ya)>270.00 ) || (MG_dir(xa,ya)<=315.00 && MG_dir(xa,ya)>270.00))  
%                 DC_quad(xa,ya)=7;
%                 MG_quad(xa,ya)=7;      
%           elseif((DC_dir(xa,ya)<=360.00 && DC_dir(xa,ya)>315.00 ) || (MG_dir(xa,ya)<=360.00 && MG_dir(xa,ya)>315.00))  
%                 DC_quad(xa,ya)=8;
%                 MG_quad(xa,ya)=8;     
%                 
%          end      
%          
%      end
%  end
         
 
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

%Unused GUI Elements
function varargout = ImageCompressionGUI(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageCompressionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageCompressionGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function pushbutton_upload_CreateFcn(hObject, eventdata, handles)
function figure1_CreateFcn(hObject, eventdata, handles)
function ImageCompressionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.UserMessages = hObject;
guidata(hObject, handles);
function varargout = ImageCompressionGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.UserMessages;
function text_stage1_Callback(hObject, eventdata, handles)
function text_stage1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function text_stage2_Callback(hObject, eventdata, handles)
function text_stage2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function text_stage3_Callback(hObject, eventdata, handles)
function text_stage3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DR_Text_Callback(hObject, eventdata, handles)
function DR_Text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function CR_Text_Callback(hObject, eventdata, handles)
function CR_Text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function CSlider_Callback(hObject, eventdata, handles)
    sliderValue = get(handles.CSlider, 'Value');  
    set(handles.SValue, 'String', num2str(round(sliderValue)))
function CSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function SValue_Callback(hObject, eventdata, handles)
function SValue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function printMsg_Callback(hObject, eventdata, handles)
function printMsg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SVDNumber_Callback(hObject, eventdata, handles)
function SVDNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SVDNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SVDSlider_Callback(hObject, eventdata, handles)
    sliderValue = get(handles.SVDSlider, 'Value');  
    set(handles.SVDNumber, 'String', num2str(round(sliderValue)))
function SVDSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SVDSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider2_Callback(hObject, eventdata, handles)
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%Main Functions
function pushbutton_upload_Callback(hObject, eventdata, handles)  
    global OrigfileSize %Button for uploading image
    [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp', 'Image Files (*.jpg, *.png, *.bmp)'}, 'Select an Image File');
    if isequal(filename,0) || isequal(pathname,0)
        return;
    else
        fullfilename = fullfile(pathname, filename);
        originalImage = imread(fullfilename);
        axes(handles.Image1); %Display uploaded image
        imshow(originalImage);
        title('Uploaded Image');
        OrigfileInfo = dir(fullfilename);
        OrigfileSize = OrigfileInfo.bytes; %Get image size and store globally
        set(handles.text_stage1, 'String', ['Original Image']);
        handles.originalImage = originalImage;
        guidata(hObject, handles);
    end

    %loads a UIgetfile window and returns if they dont select a file.
    %If they select a file then it gets stored in "originalImage"
    %Also the bytes is calculated using .bytes command and stored globally
    %to be refrenced in other functions

function pushbutton_compress_Callback(hObject, eventdata, handles) 
    global OrigfileSize
    
    if isfield(handles, 'originalImage')
        originalImage = handles.originalImage;

        %Get slider value for compression amount
        compressionSliderValue = get(handles.CSlider, 'Value');
        CompressionAmount = round(compressionSliderValue);

        outputFileName = 'compressed_image.jpg';
        outputFilePath = fullfile(pwd, outputFileName);

        tic
        %Compress and save the image
        imwrite(originalImage, outputFilePath, 'Quality', CompressionAmount);
        toc
        
        %Get file info for compressed image
        CompfileInfo = dir(outputFilePath);
        if isempty(CompfileInfo)
            msgbox('Compressed image not found.', 'Error', 'error');
            return;
        end
        CompfileSize = CompfileInfo.bytes;

        %Display compressed image and size
        CompressedImage = imread(outputFilePath);
        axes(handles.Image2);
        imshow(CompressedImage);
        title('Compressed Image');
        set(handles.text_stage2, 'String', ['Compressed Image']);

        %Store compressed image and update compression ratio
        handles.CompressedImage = CompressedImage;
        guidata(hObject, handles);

        mse = immse(originalImage, CompressedImage);

        updateCompressionRatio(OrigfileSize, CompfileSize, mse, handles)
    else
        msgbox('Please upload an image first.', 'Error', 'error');
    end

    %Function uses inbuilt compression in matlab using "imwrite" and passes
    %in the "Qaulity" parameter. This compresses the file along with a
    %value which is calculated using the slider. 

function updateCompressionRatio(originalSize, compressedSize, mse, handles)
    % ompression ratio calculation
    format short;
    compressionRatio = originalSize / compressedSize;

    displayString = sprintf('Original Size: %d\nCompressed Size: %d\nCompression Ratio: %.2f\nMSE: %.2f', ...
    originalSize, compressedSize, compressionRatio, mse);
   
    %Update Info Box
    set(handles.InfoBoxText, 'String', displayString);

    %update compression ratio just displays all the information in a box at
    %the end. including original size, new size, the ratio between them and
    %an mse value if needed for comparison (lossy)

function svdcompress_Callback(hObject, ~, handles)
    global OrigfileSize
    if isfield(handles, 'originalImage')
        originalImage = handles.originalImage;

        outputFileName = 'SVDcompressed_image.jpg';
        outputFilePath = fullfile(pwd, outputFileName);

        tic

        h = waitbar(0, 'Performing SVD Compression...');

        %Separate colour channels
        R = originalImage(:, :, 1);
        G = originalImage(:, :, 2);
        B = originalImage(:, :, 3);

        k = round(get(handles.SVDSlider, 'Value'));

        %Perform SVD for each channel
        [U_R, S_R, V_R] = svd(double(R));
        [U_G, S_G, V_G] = svd(double(G));
        [U_B, S_B, V_B] = svd(double(B));

        %Keep all singular values 
        k_R = size(S_R, 1);
        k_G = size(S_G, 1);
        k_B = size(S_B, 1);

        %Compression process for each colour channel
        compressed_R = uint8(U_R(:, 1:k) * S_R(1:k, 1:k) * V_R(:, 1:k)');
        waitbar(0.33, h, 'Performing SVD Compression - Red Channel...');
        compressed_G = uint8(U_G(:, 1:k) * S_G(1:k, 1:k) * V_G(:, 1:k)');
        waitbar(0.66, h, 'Performing SVD Compression - Green Channel...');
        compressed_B = uint8(U_B(:, 1:k) * S_B(1:k, 1:k) * V_B(:, 1:k)');
        waitbar(1, h, 'Performing SVD Compression - Blue Channel...');

        %Combine the compressed colour channels
        compressedImage = cat(3, compressed_R, compressed_G, compressed_B);

        %Save the compressed image
        imwrite(compressedImage, outputFilePath);

        %Get the compressed image file size
        CompfileInfo = dir(outputFilePath);
        CompfileSize = CompfileInfo.bytes;

        close(h);

        %Display the compressed image and size
        axes(handles.Image2);
        imshow(compressedImage);
        title('SVD Compressed Image');
        set(handles.text_stage2, 'String', ['SVD Compressed Image']);

        %Store compressed image and update compression ratio
        handles.CompressedImage = compressedImage;
        guidata(hObject, handles);

        mse = immse(originalImage, compressedImage);

        updateCompressionRatio(OrigfileSize, CompfileSize, mse, handles)

        toc
    else
        msgbox('Please upload an image first.', 'Error', 'error');
    end

    %svd splits the image down into three matricies
    %the first is how the rows relate to each other
    %the second is how the columns relate to each other
    %the third is the weight of each to show whats more important

    %for RGB we did it for each channel. "K" value only keeps the most
    %important values. More you keep the less information is lost but the
    %lower the compression rate.


function RLC_Callback(hObject, eventdata, handles)
    if isfield(handles, 'originalImage')
        set(handles.printMsg, 'String', 'Performing RLE Compression...');
        originalImage = handles.originalImage;
        
        %RUns RLE compression
        rleCompression(originalImage);

        % Display status in gui box
        set(handles.printMsg, 'String', 'RLE Compression Finished');
    else
        msgbox('Please load an image first.', 'Error', 'error');
    end

    %button function to just call the rlecompression method

function rleCompression(image)
    tic
    %Convert image to 1D array
    input_data = image(:)';
    
    %Initialise variables
    n = length(input_data);
    encoded_image = zeros(1, 2 * n);
    
    %Loop through input data
    count = 1;
    for i = 2:n
        if input_data(i) == input_data(i-1) %numbers are the same
            count = count + 1; %so add 1 to count
        else 
            encoded_image(2*i-3:2*i-2) = [input_data(i-1), count];  %updates the encoded array to the final result of that count
            count = 1; %resets the count
        end
        
        % Display progress of SVD
        if mod(i, 10000) == 0
            fprintf('Encoded %d elements, %d elements remaining\n', i, n - i);
        end
    end
    
    %Handle last run
    encoded_image(2*n-1:2*n) = [input_data(end), count];
    
    %Trim excess zeros
    encoded_image = encoded_image(encoded_image ~= 0);
    
    %Save encoded data
    save('encoded_data.mat', 'encoded_image');

    toc

function RLD_Callback(hObject, eventdata, handles)
    global OrigfileSize
    if isfield(handles, 'originalImage')
        if exist('encoded_data.mat', 'file') == 2
            set(handles.printMsg, 'String', 'Performing RLE Decompression...');
    
            %Load original image
            originalImage = handles.originalImage;
            
            %Run decompression function
            decoded_image = rleDecompression(originalImage);
    
            %Display decompressed image
            axes(handles.Image2);
            imshow(decoded_image, []); 
            title('Decompressed Image');
    
            handles.decompressedImage = decoded_image;
            guidata(hObject, handles);
    
            set(handles.printMsg, 'String', 'RLE Decompression Finished');
            fileInfo = dir('encoded_data.mat');
            compressedSize = fileInfo.bytes;

            set(handles.text_stage2, 'String', ['RLE decompressed Image']);

            updateCompressionRatio(OrigfileSize, compressedSize, 0, handles)
        else
            set(handles.printMsg, 'String', 'No compressed file found');
            return;
        end
    else
        msgbox('Please load an image first.', 'Error', 'error');
    end

    %rle decompression button to call the rleDecompression method
    %also updates the text boxes and calls the calculation method

function decoded_image = rleDecompression(originalImage)
    %Load matlab array data
    load("encoded_data.mat");
    %Initialise variables
    decoded_data = [];
    
    %Loop through the encoded data
    for i = 1:2:length(encoded_image)
        value = encoded_image(i);
        count = encoded_image(i+1);
        
        % Append the decoded data (repeated values)
        decoded_data = [decoded_data, repmat(value, 1, count)];
    end
    
    %Expected elements in original
    expected_num_elements = numel(originalImage);
    
    %Decoded_data matches the original image
    decoded_data = [decoded_data, zeros(1, expected_num_elements - length(decoded_data))];
    
    %Reshape array to the original size
    decoded_image = reshape(decoded_data, size(originalImage));

    %Save the image 
    imwrite(decoded_image, 'decoded_image.bmp'); 

    %inverses the RLE compression back into an image matrix. saves the
    %image to be returned in the button function.

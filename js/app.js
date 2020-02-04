//webkitURL is deprecated but nevertheless
URL = window.URL || window.webkitURL;

var gumStream; 						//stream from getUserMedia()
var rec; 							//Recorder.js object
var input; 							//MediaStreamAudioSourceNode we'll be recording

// shim for AudioContext when it's not avb. 
var AudioContext = window.AudioContext || window.webkitAudioContext;
var audioContext //audio context to help us record

var recordButton = document.getElementById("recordButton");
var stopButton = document.getElementById("stopButton");
var pauseButton = document.getElementById("pauseButton");

//add events to those 2 buttons
recordButton.addEventListener("click", startRecording);
stopButton.addEventListener("click", stopRecording);
pauseButton.addEventListener("click", pauseRecording);

function startRecording() {
	console.log("recordButton clicked");

	/*
		Simple constraints object, for more advanced audio features see
		https://addpipe.com/blog/audio-constraints-getusermedia/
	*/
    
    var constraints = { audio: true, video:false }

 	/*
    	Disable the record button until we get a success or fail from getUserMedia() 
	*/

	recordButton.disabled = true;
	stopButton.disabled = false;
	pauseButton.disabled = false

	/*
    	We're using the standard promise based getUserMedia() 
    	https://developer.mozilla.org/en-US/docs/Web/API/MediaDevices/getUserMedia
	*/

	navigator.mediaDevices.getUserMedia(constraints).then(function(stream) {
		console.log("getUserMedia() success, stream created, initializing Recorder.js ...");

		/*
			create an audio context after getUserMedia is called
			sampleRate might change after getUserMedia is called, like it does on macOS when recording through AirPods
			the sampleRate defaults to the one set in your OS for your playback device

		*/
		audioContext = new AudioContext({sampleRate: 44100});

		//update the format 
		document.getElementById("formats").innerHTML="Format: 1 channel pcm @ "+audioContext.sampleRate/1000+"kHz"

		/*  assign to gumStream for later use  */
		gumStream = stream;
		
		/* use the stream */
		input = audioContext.createMediaStreamSource(stream);

		/* 
			Create the Recorder object and configure to record mono sound (1 channel)
			Recording 2 channels  will double the file size
		*/
		rec = new Recorder(input,{numChannels:1, callback: testCallback})

		//start the recording process
		//rec.getBuffer(100);

		//var xx = 2;

		rec.record()

		console.log("Recording started");

	}).catch(function(err) {
	  	//enable the record button if getUserMedia() fails
    	recordButton.disabled = false;
    	stopButton.disabled = true;
    	pauseButton.disabled = true
	});
}

function testCallback(blob){
	console.log("callback");
	var newSource = audioContext.createBufferSource();
	var newBuffer = audioContext.createBuffer( 1, blob.size, audioContext.sampleRate );
	newBuffer.getChannelData(0).set(blob);
	//var newBuffer = audioContext.createBuffer( 2, buffers[0].length, audioContext.sampleRate );
	//newBuffer.getChannelData(0).set(buffers[0]);
	//newBuffer.getChannelData(1).set(buffers[1]);
	newSource.buffer = newBuffer;

	newSource.connect( audioContext.destination );
	newSource.start(0);
	console.log("end callback");
}

function pauseRecording(){
	console.log("pauseButton clicked rec.recording=",rec.recording );
	if (rec.recording){
		//pause
		//var buff = rec.getBuffer([testCallback]);
		rec.exportWAV();
		//console.log(buff);
		rec.stop();
		pauseButton.innerHTML="Resume";
	}else{
		//resume
		rec.record()
		pauseButton.innerHTML="Pause";

	}
}

function stopRecording() {
	console.log("stopButton clicked");

	//disable the stop button, enable the record too allow for new recordings
	stopButton.disabled = true;
	recordButton.disabled = false;
	pauseButton.disabled = true;

	//reset button just in case the recording is stopped while paused
	pauseButton.innerHTML="Pause";
	
	//tell the recorder to stop the recording
	rec.stop();

	rec.getBuffer(processBuffer);

	//stop microphone access
	gumStream.getAudioTracks()[0].stop();

	//create the wav blob and pass it on to createDownloadLink
	rec.exportWAV(createDownloadLink);

	// Get data byte size, allocate memory on Emscripten heap, and get pointer



}

function createDownloadLink(blob) {
	var newSource = audioContext.createBufferSource();
	var newBuffer = audioContext.createBuffer( 1, blob.size, audioContext.sampleRate );
	newBuffer.getChannelData(0).set(blob);
	newSource.buffer = newBuffer;
	
	var url = URL.createObjectURL(blob);
	var au = document.createElement('audio');
	var li = document.createElement('li');
	var link = document.createElement('a');

	//name of .wav file to use during upload and download (without extendion)
	var filename = new Date().toISOString();

	//add controls to the <audio> element
	au.controls = true;
	au.src = url;

	//save to disk link
	link.href = url;
	link.download = filename+".wav"; //download forces the browser to donwload the file using the  filename
	link.innerHTML = "Save to disk";

	//add the new audio element to li
	li.appendChild(au);
	
	//add the filename to the li
	li.appendChild(document.createTextNode(filename+".wav "))

	//add the save to disk link to li
	li.appendChild(link);
	
	//upload link
	var upload = document.createElement('a');
	upload.href="#";
	upload.innerHTML = "Upload";
	upload.addEventListener("click", function(event){
		  var xhr=new XMLHttpRequest();
		  xhr.onload=function(e) {
		      if(this.readyState === 4) {
		          console.log("Server returned: ",e.target.responseText);
		      }
		  };
		  var fd=new FormData();
		  fd.append("audio_data",blob, filename);
		  xhr.open("POST","upload.php",true);
		  xhr.send(fd);
	})
	li.appendChild(document.createTextNode (" "))//add a space in between
	li.appendChild(upload)//add the upload link to li

	//add the li element to the ol
	recordingsList.appendChild(li);
}


function processBuffer(buff) {
	/** 2D array buff[channel][index] **/
	console.log(buff);
	var audiobuff = buff[0]; // 1 channel only
	console.log(audiobuff.length);
	var top_return = 1024; // the number of samples that will be wasted to make room for the return buffer at the top of the audio array

	// pass float array, code copied from:
	// https://bl.ocks.org/jonathanlurie/e4aaa37e2d9c317ce44eae5f6011495d
	// Get data byte size, allocate memory on Emscripten heap, and get pointer

	// Import function from Emscripten generated file
    processAudioData = Module.cwrap('processAudioData', 'number', ['number', 'number','number'] );
    // function return fftargmax, float * ptr to wav array, wav length, array-return-length


    var nDataBytes = (audiobuff.length * audiobuff.BYTES_PER_ELEMENT); // 4 bytes per element for float
    var dataPtr = Module._malloc(nDataBytes); // this is emcripten specific, _malloc alloactes inside the em heap

    // Copy data to Emscripten heap (directly accessed from Module.HEAPU8)
    var dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
    dataHeap.set(new Uint8Array(audiobuff.buffer)); 
    // call function and get pitch
    //console.log(buff.length);
    //var func_in = new Uint8Array(dataHeap.buffer, dataHeap.byteOffset, audiobuff.length);
    var numblocks = processAudioData(dataHeap.byteOffset, audiobuff.length,top_return); // call emscripten function
    //var result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, audiobuff.length); // copy buffer and convert to float
    var result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset+4*(audiobuff.length-top_return), top_return); // copy top of buffer and convert to float
	//var func_out = [...result.slice(result.length-top_return)]; // top of buffer is return values

    console.log(" numblocks = "); //
    console.log(numblocks);

  // Free memory
    Module._free(dataHeap.byteOffset);
	var X = makeArr(1,numblocks,numblocks);

//*********** PLOT *************/
var trace1 = {
  x: X,
  y: result.slice(1,numblocks),
  //type: 'scatter'
  type: 'line'

};

// var trace2 = {
//   x: X,
//   y: [16, 5, 11, 9],
//   type: 'scatter'
// };
//var data = [trace1];

Plotly.newPlot('myDiv', [trace1]);

	var xx = 10;


}


function makeArr(startValue, stopValue, numpoints) {
  var arr = [];
  var step = (stopValue - startValue) / (numpoints - 1);
  for (var i = 0; i < numpoints; i++) {
    arr.push(startValue + (step * i));
  }
  return arr;
}


/** vamshi demo code
function multiAsync() {
	loadFile('test.wav', (err, data) => {
		if(err) {
			console.log('Error loading file');
			return;
		}

		a ={
			foo: 'bar',
			abc: 123,
			doSomething: () => console.log('hello')
		};

		
		{
			array1: [...],
			status: 'ok'
		}
		

	})
}
**/

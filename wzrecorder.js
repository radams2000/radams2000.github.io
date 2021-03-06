function WzRecorder(config) {
// https://www.gmass.co/blog/record-audio-mobile-web-page-ios-android/#comment-306880
    config = config || {};

    var self = this;
    var audioInput;
    var audioNode;
    var bufferSize = config.bufferSize || 4096;
    var recordedData = [];
    var recording = false;
    var recordingLength = 0;
	var startDate;
	var audioCtx;
    processAudioData = Module.cwrap('processAudioData', 'number', ['number', 'number','number','number'] );
    var nDataBytes = (4096 * 4);
    var dataPtr = Module._malloc(nDataBytes); // this is emcripten specific, _malloc alloactes inside the em heap
    var dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);


    var plotArrayLength = 60; // about 3 seconds of past history??
    var plotArray = [];
    var plotX = [];

    var PLOTYLOW = 0;
    var PLOTYHI = 1500;

    for(var i = 0; i < plotArrayLength; i++) {
        plotArray[i] = 50; // low limit for the log plot
        plotX[i] = i;
    }
    // plotArray[0] = 1500; // hi limit for log plot
    // plotArray[1] = 1500;
    var logLo = Math.log10(50);
    var logHi = Math.log10(3000);

    var layout = {
        xaxis: {range: [1, 60]},
        yaxis: {type: 'log', range: [logLo, logHi]}

        };

    

        //yaxis: {range: [0, PLOTYHI]},

    var data = { // plot all 0's to start; react is the same as newplot but faster
    x: plotX,
    y: plotArray,
    mode: 'lines',
    line: {color: '#80CAF6'}
    };

    Plotly.newPlot('myDiv', [data],layout);

 // Plotly.restyle('myDiv', 
 //     layout
 // );

    // var top_return = 2; // the number of samples that will be wasted to make room for the return buffer at the top of the audio array
    // // pass float array, code copied from:
    // // https://bl.ocks.org/jonathanlurie/e4aaa37e2d9c317ce44eae5f6011495d
    // // Get data byte size, allocate memory on Emscripten heap, and get pointer
    // // Import function from Emscripten generated file
    // processAudioData = Module.cwrap('processAudioData', 'number', ['number', 'number','number'] );
    // // function return, float * ptr to wav array, wav length, array-return-length
    // var nDataBytes = (4096 * 4); // 4 bytes per element for float
    // var dataPtr = Module._malloc(nDataBytes); // this is emcripten specific, _malloc alloactes inside the em heap
    // var dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);

    var audioprocess_counter_GL = 0;

    //console.log("here");
	
	this.toggleRecording = function()
	{
		recording ? self.stop() : self.start();
	}
	

    this.start = function() {
        document.getElementById('record').innerText = "Running, click to stop";
        var btn1 = document.getElementById('record');
        btn1.style.color = "red";
        //document.body.style.background = 'red'; 
		// reset any previous data
		recordedData = [];
		recordingLength = 0;
		
		// webkit audio context shim
		audioCtx = new (window.AudioContext || window.webkitAudioContext)();

        if (audioCtx.createJavaScriptNode) {
            audioNode = audioCtx.createJavaScriptNode(bufferSize, 1, 1);
        } else if (audioCtx.createScriptProcessor) {
            audioNode = audioCtx.createScriptProcessor(bufferSize, 1, 1);
        } else {
            throw 'WebAudio not supported!';
        }

        audioNode.connect(audioCtx.destination);

        navigator.mediaDevices.getUserMedia({audio: true})
            .then(onMicrophoneCaptured)
            .catch(onMicrophoneError);
    };

    this.stop = function() {
        // pass a function to stopRecording; the function is defined inside the call to stopRecording
        document.getElementById('record').innerText = "Paused, click to start";
        var btn2 = document.getElementById('record');
        btn2.style.color = "black";
        stopRecording(function(blob) {
			self.blob = blob;
			config.onRecordingStop && config.onRecordingStop(blob);
        });
    };
	
	this.upload = function (url, params, callback) { // not used
        var formData = new FormData();
        formData.append("audio", self.blob, config.filename || 'recording.wav');
        
        for (var i in params)
            formData.append(i, params[i]);

        var request = new XMLHttpRequest();
        request.upload.addEventListener("progress", function (e) {
            callback('progress', e, request);
        });
        request.upload.addEventListener("load", function (e) {
            callback('load', e, request);
        });

		request.onreadystatechange = function (e) {
			var status = 'loading';
			if (request.readyState == 4)
			{
				status = request.status == 200 ? 'done' : 'error';
			}
			callback(status, e, request);
		};
  
        request.open("POST", url);
        request.send(formData);
    };


    function stopRecording(callback) {
        // stop recording
        recording = false;

        // to make sure onaudioprocess stops firing
		window.localStream.getTracks().forEach( (track) => { track.stop(); });
        audioInput.disconnect();
        audioNode.disconnect();

        //console.log(recordedData.length);
        // total recording length in samples is recordingLength, set in onAudioProcess 
        // note recordedData.length gives the number of 4K sample blocks, whereas
        // recordedData(n).length gives the length of the nth sample block (4K usually)
        // recordedData is a collection of 4K-sample blocks
        var buff = mergeBuffers(recordedData,recordingLength); // convert to a single large array
        //console.log(buff.length);
        //console.log(buff);
        

        //processBuffer(buff);
		
        exportWav({
            sampleRate: sampleRate,
            recordingLength: recordingLength,
            data: recordedData
        }, function(buffer, view) {
            self.blob = new Blob([view], { type: 'audio/wav' });
            callback && callback(self.blob);
        });
    }


    function onMicrophoneCaptured(microphone) {

		if (config.visualizer)
			visualize(microphone);
		
		// save the stream so we can disconnect it when we're done
		window.localStream = microphone;

        audioInput = audioCtx.createMediaStreamSource(microphone);
        audioInput.connect(audioNode);

        audioNode.onaudioprocess = onAudioProcess;

        recording = true;
		self.startDate = new Date();
		
		config.onRecordingStart && config.onRecordingStart();
        //audioContxtOptions.sampleRate = 44100;
		sampleRate = audioCtx.sampleRate;
    }

    function onMicrophoneError(e) {
		console.log(e);
		alert('Unable to access the microphone.');
    }

    function onAudioProcess(e) {
        if (!recording) {
            return;
        }
        var num_results = 2; // C code takes in 4K data and produces 1 results (we leave some extra in case)
        audioprocess_counter_GL++;

        // don't accumulate a huge buffer, we just need a block at a time
        //recordedData.push(new Float32Array(e.inputBuffer.getChannelData(0)));
        //recordingLength += bufferSize; // 4K buffer size


        //if (audioprocess_counter_GL == 10) {
            var xxx = 0; // easy debug place
            //var audiobuff2 = e.inputBuffer.getChannelData(0);
            //var audiobuff2 = e.inputBuffer;
            //console.log(audiobuff2);

            var mydata = new Float32Array(e.inputBuffer.getChannelData(0));

            //var top_return = 2; // the number of samples that will be wasted to make room for the return buffer at the bottom of the audio array
            // pass float array, code copied from:
            // https://bl.ocks.org/jonathanlurie/e4aaa37e2d9c317ce44eae5f6011495d
            // Get data byte size, allocate memory on Emscripten heap, and get pointer
            // Import function from Emscripten generated file
            // processAudioData = Module.cwrap('processAudioData', 'number', ['number', 'number','number'] );
            // var nDataBytes = (4096 * 4), because there are 4 bytes per single-P float;
            // var dataPtr = Module._malloc(nDataBytes); // this is emcripten specific, _malloc alloactes inside the em heap
            // var dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
            
            dataHeap.set(new Uint8Array(mydata.buffer)); 
            // call emscripten function
            var numblocks = processAudioData(dataHeap.byteOffset, bufferSize,sampleRate); 
            // old var numblocks = processAudioData(dataHeap.byteOffset, bufferSize,top_return,sampleRate); 
            // copy bottom 2 floats from buffer to get result (only 0 used for now)
            var result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, num_results); 
            Module._free(dataHeap.byteOffset);
            //console.log(result);
        //}

        
        // if (audioprocess_counter_GL == 20) {
        //     console.log(numblocks);
        // }

        self.recordingLength = recordingLength;
		self.duration = new Date().getTime() - self.startDate.getTime();

		config.onRecording && config.onRecording(self.duration);

        // update the streaming plot
        //if( (result[0] > 0.0) && (result[3] > 0.0)) { // 0 result is a "flag" to stop the plot from scrolling
        if(result[0] > 0.0)  { // A 0 result is a "flag" to stop the plot from scrolling

            plotArray = plotArray.concat(result[0]);
            plotArray.splice(0, 1);

           // plotArray = plotArray.concat(result[3]);
            // plotArray.splice(0, 1);
            Plotly.update('myDiv', {y: [plotArray]},layout);

            // Plotly.relayout('myDiv', {
            //     //'yaxis.range':[0,PLOTYHI],
            //     'xaxis.autorange': false,
            //     'yaxis.autorange': false
            //     });


        }

    }




    function mergeBuffers(recBuffers, recLength) {
        var result = new Float32Array(recLength);
        var offset = 0;
        for (var i = 0; i < recBuffers.length; i++) {
            result.set(recBuffers[i], offset);
            offset += recBuffers[i].length;
        }
        return result;
    }



    function makeArr(startValue, stopValue, numpoints) {
      var arr = [];
      var step = (stopValue - startValue) / (numpoints - 1);
      for (var i = 0; i < numpoints; i++) {
        arr.push(startValue + (step * i));
      }
      return arr;
    }



	
	function visualize(stream) {
		var canvas = config.visualizer.element;
		if (!canvas)
			return;
			
		var canvasCtx = canvas.getContext("2d");
		var source = audioCtx.createMediaStreamSource(stream);

		var analyser = audioCtx.createAnalyser();
		analyser.fftSize = 2048;
		var bufferLength = analyser.frequencyBinCount;
		var dataArray = new Uint8Array(bufferLength);

		source.connect(analyser);

		function draw() {
			// get the canvas dimensions
			var width = canvas.width, height = canvas.height;

			// ask the browser to schedule a redraw before the next repaint
			requestAnimationFrame(draw);

			// clear the canvas
			canvasCtx.fillStyle = config.visualizer.backcolor || '#fff';
			canvasCtx.fillRect(0, 0, width, height);

			if (!recording)
				return;
			
			canvasCtx.lineWidth = config.visualizer.linewidth || 2;
			canvasCtx.strokeStyle = config.visualizer.forecolor || '#f00';

			canvasCtx.beginPath();

			var sliceWidth = width * 1.0 / bufferLength;
			var x = 0;

			
			analyser.getByteTimeDomainData(dataArray);

			for (var i = 0; i < bufferLength; i++) {
			
				var v = dataArray[i] / 128.0;
				var y = v * height / 2;

				i == 0 ? canvasCtx.moveTo(x, y) : canvasCtx.lineTo(x, y);
				x += sliceWidth;
			}
		
			canvasCtx.lineTo(canvas.width, canvas.height/2);
			canvasCtx.stroke();
		}
		
		draw();
	}
	
    function exportWav(config, callback) {
        function inlineWebWorker(config, cb) {

            var data = config.data.slice(0);
            var sampleRate = config.sampleRate;          
			data = joinBuffers(data, config.recordingLength);
		
            function joinBuffers(channelBuffer, count) {
                var result = new Float64Array(count);
                var offset = 0;
                var lng = channelBuffer.length;

                for (var i = 0; i < lng; i++) {
                    var buffer = channelBuffer[i];
                    result.set(buffer, offset);
                    offset += buffer.length;
                }

                return result;
            }

            function writeUTFBytes(view, offset, string) {
                var lng = string.length;
                for (var i = 0; i < lng; i++) {
                    view.setUint8(offset + i, string.charCodeAt(i));
                }
            }

            var dataLength = data.length;

            // create wav file
            var buffer = new ArrayBuffer(44 + dataLength * 2);
            var view = new DataView(buffer);
			
            writeUTFBytes(view, 0, 'RIFF'); // RIFF chunk descriptor/identifier
            view.setUint32(4, 44 + dataLength * 2, true); // RIFF chunk length
            writeUTFBytes(view, 8, 'WAVE'); // RIFF type
            writeUTFBytes(view, 12, 'fmt '); // format chunk identifier, FMT sub-chunk
            view.setUint32(16, 16, true); // format chunk length
            view.setUint16(20, 1, true); // sample format (raw)
            view.setUint16(22, 1, true); // mono (1 channel)
            view.setUint32(24, sampleRate, true); // sample rate
            view.setUint32(28, sampleRate * 2, true); // byte rate (sample rate * block align)
            view.setUint16(32, 2, true); // block align (channel count * bytes per sample)
            view.setUint16(34, 16, true); // bits per sample
            writeUTFBytes(view, 36, 'data'); // data sub-chunk identifier
            view.setUint32(40, dataLength * 2, true); // data chunk length

            // write the PCM samples
            var index = 44;
            for (var i = 0; i < dataLength; i++) {
                view.setInt16(index, data[i] * 0x7FFF, true);
                index += 2;
            }

            if (cb) {
                return cb({
                    buffer: buffer,
                    view: view
                });
            }

            postMessage({
                buffer: buffer,
                view: view
            });
        }

        var webWorker = processInWebWorker(inlineWebWorker);

        webWorker.onmessage = function(event) {
            callback(event.data.buffer, event.data.view);

            // release memory
            URL.revokeObjectURL(webWorker.workerURL);
        };

        webWorker.postMessage(config);
    }

    function processInWebWorker(_function) {
        var workerURL = URL.createObjectURL(new Blob([_function.toString(),
            ';this.onmessage = function (e) {' + _function.name + '(e.data);}'
        ], {
            type: 'application/javascript'
        }));

        var worker = new Worker(workerURL);
        worker.workerURL = workerURL;
        return worker;
    }
}
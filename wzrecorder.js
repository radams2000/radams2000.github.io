function WzRecorder(config) {

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
	
	this.toggleRecording = function()
	{
		recording ? self.stop() : self.start();
	}
	

    this.start = function() {

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

        console.log(recordedData.length);
        // total recording length in samples is recordingLength, set in onAudioProcess 
        // note recordedData.length gives the number of 4K sample blocks, whereas
        // recordedData(n).length gives the length of the nth sample block (4K usually)
        // recordedData is a collection of 4K-sample blocks
        var buff = mergeBuffers(recordedData,recordingLength); // convert to a single large array
        //console.log(buff.length);
        //console.log(buff);
        processBuffer(buff);
		
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

        recordedData.push(new Float32Array(e.inputBuffer.getChannelData(0)));
        recordingLength += bufferSize;

        self.recordingLength = recordingLength;
		self.duration = new Date().getTime() - self.startDate.getTime();

		config.onRecording && config.onRecording(self.duration);
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


    function processBuffer(audiobuff) {
        /** 2D array buff[channel][index] **/
        //console.log(buff);
        //var audiobuff = buff[0]; // 1 channel only
        //console.log(audiobuff.length);
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
        console.log(" C-code results in top 1K of heap  = "); 
        console.log(result);

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
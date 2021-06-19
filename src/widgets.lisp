(in-package :deliscope)

(defclass parsers ()
  ((parsers :initform (make-hash-table :test 'eq :thread-safe t)
            :initarg :parsers
            :accessor parsers))
  )

(defun estimate-sequences-in-file (file-name &key num verbose)
  (declare (ignore verbose))
  (unless (probe-file file-name)
    (error "Could not find file: ~a" file-name))
  (if num
      num
      (multiple-value-bind (num-lines file-pos file-size)
          (core:count-lines-in-file file-name 1000000)
          (let ((estimate-lines (floor (* (/ file-size file-pos) num-lines))))
            (floor (/ estimate-lines 4))))))

(defparameter *data-directory* nil)
(defun set-data-directory (directory)
  (setf directory (uiop:ensure-directory-pathname directory))
  (unless (probe-file directory)
    (error "The directory ~a does not exist" directory)
    (setf *data-directory* directory)))

(defun create-sequence-parser (name files &key num output-file-name (overwrite t))
  (unless output-file-name
    (setf output-file-name (make-pathname :name (format nil "~a-results" (string-downcase (string name)))
                                          :type "dat")))
  (unless overwrite
    (when (probe-file output-file-name)
      (error "Saving results will fail because the file ~a already exists" output-file-name)))
  (let ((absolute-files (mapcar (lambda (file)
                                  (if *data-directory*
                                      (merge-pathname file :default *data-directory*)
                                      file))
                                files)))
    (format t "Scanning the files to count the number of sequences.~%")
    (format t "This may take a few minutes - hit ii to interrupt.~%")
    (finish-output)
    (let* ((num-seq-per-file (lparallel:pmapcar (lambda (file)
                                                  (estimate-sequences-in-file file :num num :verbose nil))
                                                files))
           (parser (make-instance 'parser :name name
                                          :files absolute-files
                                          :num-sequences-per-file num-seq-per-file
                                          :pass-file-name pass-file-name
                                          :fail-file-name fail-file-name
                                          :output-file-name output-file-name
                                          :overwrite overwrite)))
      (format t "Number of sequences: ~a~%" (apply '+ num-seq-per-file))
      (finish-output)
      parser)))


(defparameter *progress* nil)
(defparameter *parser* nil)

(defclass progress ()
  ((progress-bar :initarg :progress-bar :reader progress-bar)
   (messages :initarg :messages :reader messages)
   (holder :initarg :holder :accessor holder)))

(defun parse-impl (serial-parallel parsers)
  (let* ((container (make-instance 'jw:v-box))
         (panels (loop for parser in parsers
                       collect (cw:make-threaded-task-page
                                container
                                (format nil "Task ~a" (name parser))
                                (lambda (instance action parser progress-callback)
                                  (declare (ignore action))
                                  (funcall serial-parallel
                                           parser
                                           :progress-callback
                                           (let ((last-val -1))
                                             (lambda (val msg &key done)
                                               (if done
                                                   (format t "Done.~%")
                                                   (when (> val last-val)
                                                     (setf last-val val)
                                                     (progn
                                                       (format t "Progress: ~a~%" val)
                                                       (format t "~a~%" msg)
                                                       (finish-output)
                                                       (funcall progress-callback :value val
                                                                                  :maximum 100)))))))
                                  t)
                                :parameter parser
                                :label "Click button to start."))))
    (j:display container)
    (loop for panel in panels
          do (cw:run-task panel))
    (values)))

(defun parse (&rest parsers)
  "Parse sequence files one at a time but parsers are run in parallel."
  (parse-impl 'serial-analyze parsers))

(defun parse-fast (&rest parsers)
  "Parse sequence files in parallel and parsers in parallel.
 This is experimental and depends on the SeqAn 2.0 library to be thread safe,
 which I am not absolutely certain it is.  If this crashes - use 'run'."
  (parse-impl 'parallel-analyze parsers))

(defun plot-counts (data title)
  (let* ((x-data (coerce (subseq (mapcar #'car data) 0) 'vector))
         (y-data (coerce (subseq (mapcar #'cdr data) 0) 'vector)))
    (dv:bar-plot x-data
                 y-data
                 :width 600
                 :x-title "index" :y-title "Log10 counts"
                 :title title)))

(defun plot-unfiltered-counts (analysis)
  (let ((data (sequence-counts (sequences analysis))))
    (plot-counts data (format nil "~a unfiltered counts" (name analysis)))))

(defun plot-filtered-counts (analysis)
  (let ((data (sequence-counts (filtered-sequences analysis))))
    (plot-counts data (format nil "~a filtered counts" (name analysis)))))


#|

Copyright (c) 2021, Christian E. Schafmeister

Deliscope is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

See file 'LICENSE' for full details.

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

|#

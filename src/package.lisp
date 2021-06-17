
(defpackage :deliscope
  (:use :common-lisp)
  (:nicknames :ds)
  (:shadow sequence count)
  (:export
   #:*minimum-phred-quality*
   #:sequences
   #:set-filter
   #:library-codes
   #:filter-seq-file
   #:load-encoding
   #:analyze
   #:parallel-analyze
   #:serial-analyze
   #:clear-analysis
   #:remove-rare-sequences
   #:join-codons
   #:filter-sequences
   #:bead-specific-counts
   #:sort-hits
   #:*analysis*
   #:save-results
   #:load-results
   #:multiple-redundant
   #:merge-hits



   #:compare
   #:create-sequence-parser
   #:filtered-library-codes
   #:summarize
   #:redundancies
   #:save-csv
   #:parse
   #:parse-fast
   ))
   


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

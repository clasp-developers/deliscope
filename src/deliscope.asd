(in-package :asdf-user)

(defsystem "deliscope"
  :description "Analyze a DEL DNA sequencing run"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:seqan)
  :serial t
  :components
  ((:file "package")
   (:file "filter")
   (:file "analyze")
   (:file "widgets")
   ))

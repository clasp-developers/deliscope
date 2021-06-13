
(format t "*load-pathname* -> ~a~%" (truename *load-pathname*))
(format t "Starting up del analysis~%")
(format t "deliscope asd path: ~s~%" (merge-pathnames #P"deliscope.asd" (truename *load-pathname*)))
(ql:quickload :delta-vega)
(asdf:load-asd (merge-pathnames #P"deliscope.asd" (truename *load-pathname*)))
(asdf:load-system :deliscope)
(use-package :ds)

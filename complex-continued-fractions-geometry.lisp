(defpackage :complex-continued-fractions-geometry
  (:nicknames :ccfg)
  (:use :cl :ol :iterate :moebius-transformations)
  (:export))

(in-package :complex-continued-fractions-geometry)

;;; TODO collect some of the import regions of the Schmidt CF
;;; expansion, and also the transforms used with it
(defparameter  maps
  (list :v1 (mt 1 i
                0 1)
        :v2 (mt 1 0
                (- i) 1)
        :v3 (mt (- 1 i) i
                (- i) (+ 1 i))
        :e1 (mt 1 0
                (- 1 i) i)
        :e2 (mt 1 (- i 1)
                0 i)
        :e3 (mt i 0
                0 1)
        :c (mt 1 (- i 1)
               (- 1 i) i)
        :s (mt 0 -1
               1 -1)
        :i (mt 1 0
               0 1)))

(defparameter regions
  ;; for V
  (list :v1 (line i 1)
        :v2 (projective-line-3points 0 i (/ (- 1 i)))
        :v3 (projective-line-3points 1 (+ 1 i) (/ (- 1 i)))
        :r (line 0 1)
        ;; for V^*
        :a1 (projective-line-3points 0 1 (/ (- 1 i)))
        :a2 (line 0 i)
        :a3 (line 1 i)
        :c* (projective-line-3points i (+ 1 i) (/ (- 1 i)))))

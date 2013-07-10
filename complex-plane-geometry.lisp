(defpackage :complex-plane-geometry
  (:nicknames :cgeom)
  (:shadow :transformation)
  (:shadowing-import-from :moebius-transformations #:line)
  (:use :clim :clim-lisp :ol :iterate
        :moebius-transformations)
  (:export))

(in-package :complex-plane-geometry)

(defparameter window-width 1000)
(defparameter window-height 800)

(define-application-frame complex-plane ()
  ((bound :initform 1000
          :accessor bound)
   (zoom :initform 1
         :accessor zoom)
   (objects :initarg :objects
            :initform nil
            :accessor objects))
  (:panes (canvas :application
                  :width window-width :height window-height
                  :incremental-redisplay t
                  :scroll-bars t)
          (int :interactor
               :width 400 :height window-height))
  (:layouts (default (horizontally () int canvas))))

(progn
  (defparameter pixels-per-one 100)

  (let ((h-center (floor window-width 2))
        (v-center (floor  window-height 2)))
    (defparameter transformation
      (make-3-point-transformation*
       ;; from
       0 0
       1 0
       0 1
       ;; to
       h-center v-center
       (+ h-center pixels-per-one) v-center
       h-center (- v-center pixels-per-one)))))

;; this gives a drawing region of [-50,50] x [-40, 40]

(defun complex-plane-interactor ()
  (run-frame-top-level (make-instance 'complex-plane)))

(defmacro with-canvas (&body body)
  `(let ((canvas (get-frame-pane *application-frame* 'canvas)))
     ,@body))

(define-presentation-type projective-line ())

(define-presentation-type line () :inherit-from 'projective-line)

(define-presentation-type circle () :inherit-from 'projective-line)

(defmethod add-object ((line line))
  ;; determine the edge points
  (mvbind (min max) (intersect-rectangle line :rmin -50 :rmax 50 :imin -40 :imax 40)
    (let ((point-1 (+ (basepoint line) (* min (direction line))))
          (point-2 (+ (basepoint line) (* max (direction line)))))
      (with-canvas
        (with-output-as-presentation (canvas line 'line
                                             :allow-sensitive-inferiors nil
                                             :single-box t)
          (draw-line* canvas
                      (realpart point-1) (imagpart point-1)
                      (realpart point-2) (imagpart point-2)
                      :transformation transformation))))))

(defmethod add-object :after (object)
  (pushnew object (objects *application-frame*)))


(defmethod add-object ((circle circle))
  (with-canvas
    (with-output-as-presentation (canvas circle 'circle
                                         :allow-sensitive-inferiors nil
                                         :single-box t)
      (draw-circle* (get-frame-pane *application-frame* 'canvas)
                    (realpart (center circle))
                    (imagpart (center circle))
                    (radius circle)
                    :filled nil
                    :transformation transformation))))

(define-complex-plane-command (com-hover :name "Hover test" :menu t)
    ((projective-line 'projective-line)
     (line 'line)
     (circle 'circle)))

(define-complex-plane-command (com-add-circle
                               :name "Add circle (center & radius)" :menu t)
    ((center 'number :prompt "Center")
     (radius 'real :prompt "Radius"))
  (add-object (make-instance 'circle :center center :radius radius)))

(define-complex-plane-command (com-add-line
                               :name "Add line (point & direction)" :menu t)
    ((basepoint 'number :prompt "Point")
     (direction 'number :prompt "Direction"))
  (add-object (make-instance 'line :basepoint basepoint :direction direction)))

(define-complex-plane-command (com-translate
                               :name "Translate" :menu t)
    ((object 'projective-line)
     (offset 'number))
  (add-object
   (translate offset object)))

(define-complex-plane-command (com-draw-farey-sets :name "Draw Farey sets" :menu t)
    ()
  (dolist (r ccfg:regions)
    (unless (symbolp r)
      (add-object  r))))

(define-complex-plane-command (com-clear-canvas :name "Clear Canvas" :menu t)
    ()
  (with-canvas (window-clear canvas)))

;;; TODO get presentations to work

;;; TODO filled circles,
;;; TODO filled half-planes,
;;; TODO filled triangles

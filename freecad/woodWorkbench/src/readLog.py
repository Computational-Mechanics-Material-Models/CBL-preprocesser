from pathlib import Path
import FreeCAD as App # type: ignore

def readLog(self):
    
    paraFile = self.form[0].setupFile.text()

    # Read parameters from file
    with open(Path(paraFile), "r") as f:
        for line in f:
            if "geoName" in line:
                geoName = line.split("=")[1].strip()
            if "radial_growth_rule" in line:
                radial_growth_rule = line.split("=")[1].strip()
            if "species" in line:
                species = line.split("=")[1].strip()
            elif "ring_width" in line:
                ring_width = line.split("=")[1].strip()
            elif "late_ratio" in line:
                late_ratio = line.split("=")[1].strip()
            elif "cellsize_early" in line:
                cellsize_early = line.split("=")[1].strip()
            elif "cellsize_late" in line:
                cellsize_late = line.split("=")[1].strip()
            elif "cellwallthickness_early" in line:
                cellwallthickness_early = line.split("=")[1].strip()
            elif "cellwallthickness_late" in line:
                cellwallthickness_late = line.split("=")[1].strip()
            elif "cell_length" in line:
                cell_length = line.split("=")[1].strip()

            elif "randomFlag" in line:
                randomFlag = line.split("=")[1].strip()
            elif "dist_type" in line:
                dist_type = line.split("=")[1].strip()
            elif "dist_params" in line:
                dist_params = line.split("=")[1].strip().replace(']', '').replace('[', '')
            elif "corr_l" in line:
                corr_l = line.split("=")[1].strip()
            elif "sampling_type" in line:
                sampling_type = line.split("=")[1].strip()

            elif "box_shape" in line:
                box_shape = line.split("=")[1].strip()
            elif "box_center" in line:
                box_center = line.split("=")[1].strip()
            elif "box_size" in line:
                box_size = line.split("=")[1].strip()
            elif "box_height" in line:
                box_height = line.split("=")[1].strip()
            elif "box_width" in line:
                box_width = line.split("=")[1].strip()
            elif "box_depth" in line:
                box_depth = line.split("=")[1].strip()
            elif "x_notch_size" in line:
                x_notch_size = line.split("=")[1].strip()
            elif "y_notch_size" in line:
                y_notch_size = line.split("=")[1].strip()

            elif "precrackFlag" in line:
                precrackFlag = line.split("=")[1].strip()
            elif "precrack_depth" in line:
                precrack_depth = line.split("=")[1].strip()
            elif "precrack_width" in line:
                precrack_width = line.split("=")[1].strip()

            elif "iter_max" in line:
                iter_max = line.split("=")[1].strip()
            elif "theta_min" in line:
                theta_min = line.split("=")[1].strip()
            elif "long_connector_ratio" in line:
                long_connector_ratio = line.split("=")[1].strip()

            elif "knotFlag" in line:
                knotFlag = line.split("=")[1].strip()
            elif "Uinf" in line:
                Uinf = line.split("=")[1].strip()
            if "a1" in line:
                a1 = line.split("=")[1].strip()
            elif "a2" in line:
                a2 = line.split("=")[1].strip()
            elif "m1" in line:
                m1 = line.split("=")[1].strip()
            elif "m2" in line:
                m2 = line.split("=")[1].strip()

            elif "boundaryFlag" in line:
                boundaryFlag = line.split("=")[1].strip()
            elif "flowFlag" in line:
                flowFlag = line.split("=")[1].strip()
            elif "mergeFlag" in line:
                mergeFlag = line.split("=")[1].strip()
            elif "rayFlag" in line:
                rayFlag = line.split("=")[1].strip()
            elif "inpType" in line:
                inpType = line.split("=")[1].strip()

            elif "geoFil" in line:
                geoFil = line.split("=")[1].strip()

    # Write parameters to input panel
        self.form[0].geoName.setText(geoName)
        self.form[0].radial_growth_rule.setCurrentText(radial_growth_rule)
        self.form[0].species.setCurrentText(species)
        self.form[0].ring_width.setText(ring_width)
        self.form[0].ring_ratio.setText(late_ratio)
        self.form[0].cellsize_early.setText(cellsize_early)
        self.form[0].cellsize_late.setText(cellsize_late)
        self.form[0].cellwallthickness_early.setText(cellwallthickness_early)
        self.form[0].cellwallthickness_late.setText(cellwallthickness_late)
        self.form[0].cell_length.setText(cell_length)

        self.form[0].randomFlag.setCurrentText(randomFlag)
        self.form[0].dist_types.setCurrentText(dist_type)
        self.form[0].dist_params.setText(dist_params)
        self.form[0].corr_l.setText(corr_l)
        self.form[0].sampling_type.setCurrentText(sampling_type)

        
        self.form[1].box_center.setText(box_center)
        if box_shape == 'cube':
            self.form[1].box_shape.setCurrentText('Cube')
            self.form[1].cube_size.setText(box_size)
        elif box_shape == 'rectangle':
            self.form[1].box_shape.setCurrentText('Rectangle')
            self.form[1].box_height.setText(box_height)
            self.form[1].box_width.setText(box_width)
            self.form[1].box_depth.setText(box_depth)
        elif box_shape == 'notchedsquare':
            self.form[1].box_shape.setCurrentText('Notched Square')
            self.form[1].notch_height.setText(box_height)
            self.form[1].notch_width.setText(box_width)
            self.form[1].notch_depth.setText(box_depth)
            self.form[1].x_indent_size.setText(x_notch_size)
            self.form[1].y_indent_size.setText(y_notch_size)
        if box_shape == 'input':
            self.form[1].box_shape.setCurrentText('Input')
            self.form[1].geoFile.setText(geoFil)
            self.form[1].geo_size.setText(box_size)
            self.form[1].geo_height.setText(box_height)

        self.form[1].precrackFlag.setCurrentText(precrackFlag)
        self.form[1].precrack_depth.setText(precrack_depth)
        self.form[1].precrack_width.setText(precrack_width)

        self.form[1].iter_max.setText(iter_max)
        self.form[1].theta_min.setText(theta_min)
        self.form[1].long_connector_ratio.setText(long_connector_ratio)
        self.form[1].knotFlag.setCurrentText(knotFlag)
        self.form[1].knot_flow.setText(Uinf)
        self.form[1].a1.setText(a1)
        self.form[1].a2.setText(a2)
        self.form[1].m1.setText(m1)
        self.form[1].m2.setText(m2)

        self.form[2].boundaryFlag.setCurrentText(boundaryFlag)
        self.form[2].flowFlag.setCurrentText(flowFlag)
        self.form[2].merge_operation.setCurrentText(mergeFlag)
        self.form[2].rayFlag.setCurrentText(rayFlag)
        self.form[2].inpType.setCurrentText(inpType)



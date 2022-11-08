package com.aerodynamics.airfoil_calculation;

import androidx.appcompat.app.AppCompatActivity;
import org.tensorflow.lite.Interpreter;

import android.content.Context;
import android.content.res.AssetFileDescriptor;
import android.graphics.Color;
import android.graphics.PorterDuff;
import android.os.Bundle;
import android.view.View;
import android.view.inputmethod.InputMethodManager;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

import com.google.android.material.switchmaterial.SwitchMaterial;
import com.google.firebase.analytics.FirebaseAnalytics;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Arrays;

public class AI_Inverse_problem extends AppCompatActivity {

	SwitchMaterial naca4digits, naca5digits;
	AssetFileDescriptor fileDescriptor;
	EditText clVisc, cdVisc, aOa;
	Button button;
	float clvisc;
	float cdvisc;
	float aoa;
	TextView result4, result5;
	Interpreter tflite;
	float[] min_data = new float[]{-8.0f, -1.06830704f, -0.09636456f};
	float[] max_data = new float[]{8.0f, 3.4748194f, 0.08257051f};
	float[] min_target = new float[]{0.0f, 0.0f, 1.0f};
	float[] max_target = new float[]{9.0f, 9.0f, 30.0f};

	float[] min_data5   = new float[]{-8.0f, -0.975260480f, 2.34188541e-3f};
	float[] max_data5   = new float[]{8.0f, 1.2971378f, 0.09098151f};
	float[] min_target5 = new float[]{10.0f, 1.0f};
	float[] max_target5 = new float[]{50.0f, 30.0f};

	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_ai_inverse_problem);

		FirebaseAnalytics mFirebaseAnalytics = FirebaseAnalytics.getInstance(this);

		naca4digits = findViewById(R.id.switch4);
		naca5digits = findViewById(R.id.switch5);

		button = findViewById(R.id.button);
		int defCol = button.getTextColors().getDefaultColor();

		clVisc = findViewById(R.id.editClVisc);
		clVisc.setHint("0.0");

		cdVisc = findViewById(R.id.editCdVisc);
		cdVisc.setHint("0.002");

		aOa = findViewById(R.id.editAOa);
		aOa.setHint("0.0");

		result4 = findViewById(R.id.result4);
		result5 = findViewById(R.id.result5);

		naca4digits.setOnCheckedChangeListener((buttonView, isChecked) -> {
			if (isChecked) {
				naca5digits.setChecked(false);

				button.setEnabled(true);
				button.getBackground().setColorFilter(null);
				button.setTextColor(defCol);
			}
		});
		naca5digits.setOnCheckedChangeListener((buttonView, isChecked) -> {
			if (isChecked) {
				naca4digits.setChecked(false);

				button.setEnabled(true);
				button.getBackground().setColorFilter(null);
				button.setTextColor(defCol);
			}
		});
		if (!naca4digits.isChecked() && !naca5digits.isChecked()) {
			button.setEnabled(false);
			button.getBackground().setColorFilter(Color.rgb(245, 245, 245), PorterDuff.Mode.MULTIPLY);
			button.setTextColor(Color.LTGRAY);
		}
	}

	public void onClick(View view) {

		try {
			tflite = new Interpreter(loadModelFile());
		}catch (Exception ex){
			ex.printStackTrace();
		}

		view = this.getCurrentFocus();
		if (view != null) {
			InputMethodManager imm = (InputMethodManager) getSystemService(Context.INPUT_METHOD_SERVICE);
			imm.hideSoftInputFromWindow(view.getWindowToken(), 0);
		}

		if (clVisc.getText().toString().isEmpty()) {
			clvisc = 0.0f;
		} else
			clvisc = Float.parseFloat(clVisc.getText().toString());

		if (cdVisc.getText().toString().isEmpty()) {
			cdvisc = 0.02f;
		} else
			cdvisc = Float.parseFloat(cdVisc.getText().toString());

		if (aOa.getText().toString().isEmpty()) {
			aoa = 0.0f;
		} else
			aoa = Float.parseFloat(aOa.getText().toString());

		if (naca4digits.isChecked()) {
			result4.setVisibility(View.VISIBLE);
			result5.setVisibility(View.GONE);

			float[][] prediction = doInference(aoa, clvisc, cdvisc);
			System.out.println(Arrays.deepToString(prediction));
			if(prediction[0][2] > 10.0) {
				result4.setText("The corresponding profile is: " + "NACA" + Math.round(prediction[0][0]) + Math.round(prediction[0][1]) + Math.round(prediction[0][2]));
			}
			else{
				result4.setText("The corresponding profile is: " + "NACA" + Math.round(prediction[0][0]) + Math.round(prediction[0][1]) + "0" + Math.round(prediction[0][2]));
			}
		} else if (naca5digits.isChecked()) {
			result4.setVisibility(View.GONE);
			result5.setVisibility(View.VISIBLE);

			float[][] prediction = doInference(aoa, clvisc, cdvisc);
			System.out.println(Arrays.deepToString(prediction));

			float[] values = new float[]{10, 20, 30, 40, 50};
			float prediction_1 = Math.round(prediction[0][0]);
			float distance = Math.abs(values[0] - prediction_1);
			int idx = 0;
			for(int c = 1; c < values.length; c++){
				float cdistance = Math.abs(values[c] - prediction_1);
				if(cdistance < distance){
					idx = c;
					distance = cdistance;
				}
			}
			prediction_1 = values[idx];

			if(prediction[0][1] > 10.0) {
				result5.setText("The corresponding profile is: " + "NACA 2" + Math.round(prediction_1) + Math.round(prediction[0][1]));
			}
			else{
				result5.setText("The corresponding profile is: " + "NACA 2" + Math.round(prediction_1) + "0" + Math.round(prediction[0][1]));
			}
		}
	}


	private MappedByteBuffer loadModelFile() throws IOException {
		if(naca4digits.isChecked()) {
			fileDescriptor = this.getAssets().openFd("trained_model.tflite");
		} else if(naca5digits.isChecked()){
			fileDescriptor = this.getAssets().openFd("trained_model_5.tflite");
		}
		FileInputStream inputStream = new FileInputStream(fileDescriptor.getFileDescriptor());
		FileChannel fileChannel = inputStream.getChannel();
		long startOffset = fileDescriptor.getStartOffset();
		long declareLength = fileDescriptor.getDeclaredLength();
		return fileChannel.map(FileChannel.MapMode.READ_ONLY, startOffset, declareLength);
	}

	private float[][] doInference(float aOa, float cl, float cd) {

		if(naca4digits.isChecked()) {
			aOa = (aOa - min_data[0]) / (max_data[0] - min_data[0]);
			cl = (cl - min_data[1]) / (max_data[1] - min_data[1]);
			cd = (cd - min_data[2]) / (max_data[2] - min_data[2]);
		}
		else if(naca5digits.isChecked()){
			aOa = (aOa - min_data5[0]) / (max_data5[0] - min_data5[0]);
			cl = (cl - min_data5[1]) / (max_data5[1] - min_data5[1]);
			cd = (cd - min_data5[2]) / (max_data5[2] - min_data5[2]);
		}

		float[] inputVal = new float[3];
		inputVal[0] = aOa;
		inputVal[1] = cl;
		inputVal[2] = cd;
		float[][] output4 = new float[1][3];
		float[][] output5 = new float[1][2];
		if(naca4digits.isChecked()) {
			tflite.run(inputVal, output4);
		}
		else if(naca5digits.isChecked()){
			tflite.run(inputVal, output5);
			}
		if (naca4digits.isChecked()) {
			output4[0][0] = output4[0][0] * (max_target[0] - min_target[0]) + min_target[0];
			output4[0][1] = output4[0][1] * (max_target[1] - min_target[1]) + min_target[1];
			output4[0][2] = output4[0][2] * (max_target[2] - min_target[2]) + min_target[2];
			return output4;
		}
		else if (naca5digits.isChecked()){
			output5[0][0] = output5[0][0] * (max_target5[0] - min_target5[0]) + min_target5[0];
			output5[0][1] = output5[0][1] * (max_target5[1] - min_target5[1]) + min_target5[1];
			return output5;
		}
		return output4;
	}
}
using System;
using UnityEngine;
using UnityEngine.UI;
using Input = UnityEngine.Input;

namespace NeroAdventure.Physics.Field
{
    [RequireComponent(typeof(RawImage))]
    public class FlowFieldTexture : MonoBehaviour
    {
        [SerializeField] private int length = 128;
        [SerializeField] private float viscosityRate;
        [SerializeField] private float diffuseRate;
        
        private RawImage rawImage;
        private FlowField flowField;
        private Texture2D texture;

        private float[,] flow;


        private Vector2 mousePos_prev;

        private void Awake()
        {
            flowField = new FlowField(length, viscosityRate, diffuseRate);
            texture = new Texture2D(length, length);
            
            flow = new float[length, length];
            mousePos_prev = Vector2.zero;
        }

        private void Start()
        {
            rawImage = GetComponent<RawImage>();
            rawImage.texture = texture;
        }

        private void Update()
        {
            Vector2 mousePos = GetPos();
            if (UnityEngine.Input.GetMouseButton(0))
            {
                Vector2Int pos = Vector2Int.RoundToInt(mousePos);
                
                flowField.AddDensity(pos, 1f);
                flowField.AddVelocity(pos, mousePos - mousePos_prev);
            }
            
            flowField.Step();
            Draw();
            mousePos_prev = mousePos;
        }
        
        private Vector2 GetPos()
        {
            Vector2 mousePos = UnityEngine.Input.mousePosition;
            Vector2 imagePos = rawImage.rectTransform.position;
            Vector2 div = new Vector2(rawImage.rectTransform.rect.width, rawImage.rectTransform.rect.height);
            Vector2 pos = (mousePos - imagePos + div / 2) / div * length;

            float x = Mathf.Clamp(pos.x, 0, length - 1);
            float y = Mathf.Clamp(pos.y, 0, length - 1);

            return new Vector2(x, y);
        }

        private void Draw()
        {
            float[,] flow = flowField.GetDensityField();
            
            for (int i = 0; i < length; i++)
            {
                for (int j = 0; j < length; j++)
                {
                    float density = flow[i, j];
                    Color c = new Color(density, density, density);
                    texture.SetPixel(j, i, c);
                }
            }
            
            texture.Apply();
        }
    }
}
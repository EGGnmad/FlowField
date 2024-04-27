using System;
using UnityEngine;

namespace NeroAdventure.Physics.Field
{
    public class FlowFieldBehavior : MonoBehaviour
    {
        [Header("Flow Field")]
        [SerializeField] private int length;
        [SerializeField] private float diffuseRate;
        [SerializeField] private float viscosityRate;

        private FlowField flowField;

        private void Awake()
        {
            flowField = new FlowField(length, viscosityRate, diffuseRate);
        }

        // Debug Size
        private void OnDrawGizmosSelected()
        {
            Gizmos.color = Color.cyan;
            Gizmos.DrawWireCube(transform.position, transform.localScale);
        }
    }
}